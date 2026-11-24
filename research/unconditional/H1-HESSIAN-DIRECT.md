# H1-HESSIAN-DIRECT — the expected Hessian at the soft k→0 mode: Savvidy pointwise-indefinite but expectation-positive, and the gap is *between* the two (a conditioning defect, not a sign)

> 2026-11-23 · Direct head-on attack on `G″(k→0)≥λ`, signed at the expected Hessian on the gauge-invariant soft mode. Decompose `E[H]=E[−D²]+E[2adF_B]+E[Φ″]`, verify each sign, and ask whether the Savvidy term `E[2adF_B]` eats the budget.
>
> **判词 (verdict, up front).** The crux resolves into a **true theorem plus a false reduction**.
> **(1) A true theorem.** The exact Schwinger–Dyson/Stein identity makes the **full** expected Hessian positive-semidefinite *unconditionally*: `E_μ[H_vv]=E_μ[f²]⪰0` for every measure, every β, no convexity. Hence `E[2adF_B]=E[f²]−E[−D²]−E[Φ″]≥−(E[−D²]+E[Φ″])` **exactly** — the Savvidy term provably **does not eat the budget**; the expected Hessian is genuinely expectation-positive even though pointwise it is indefinite (bottom `−g|B|`). **(2) A correction.** The claimed identity is **not** `E[−D²]=E[f²]`: `−D²` is only the gauge-covariant kinetic (Laplacian) piece; the Stein identity carries `f²` for the **full** `−D²+2adF_B+Φ″`. Mis-attributing it to `−D²` silently deletes exactly the sign-carrying (`2adF_B`) and budget (`Φ″`) pieces. **(3) A false reduction.** `G″≥λ` does **not** reduce to a single sign of `E[2adF_B]`. `G″(0)=−(log p₀)″(0)` is the **conditioned** coarse-marginal curvature; unconditionally at the `P`-symmetric vacuum `E[H]=E[f²]=Var(f)` **exactly** (equipartition), so the flat object is identically `0` and the entire gap is the **conditioning correction** (`W_geo+X`). The measured `D/S≈2%` shows the coarse marginal is `98%` Gaussian (`G″≈0.98/Var(a)>0`) — expectation-positive curvature with a `2%` margin — but `D<S` is *measured, not proven*; its proof is the conditioned force-variance bound `Var_μ(f|a{=}0)≤E[H|a{=}0]+W−λ`, which **is H1**. The gap lives **between the unconditioned expectation (provably `≥0`) and the conditioned curvature (H1)** — not between "pointwise indefinite" and "expectation-positive", both of which are settled and on the wrong side of the gap.

---

## §1 · The decomposition on the gauge-invariant soft mode — exact objects and their signs

**Setup.** Lattice SU(N) Wilson, `μ=Z⁻¹e^{−S_W}dU`, `S_W=(1/2g²)Σ_□tr F²`. Flat chart `X_ℓ=log U_ℓ∈su(N)`, `dU=J dX`, `Φ=−log J`, so `μ∝e^{−V}dX`, `V=S_W+Φ`. The **gauge-invariant soft mode** is the coarse block-uniform / zero-momentum transverse direction `v` (project `Q⊥` onto the gauge orbit; `a=⟨v,X⟩`, `|v|=1`, `v` transverse). On this mode the microscopic second variation is (PROOF-SAVVIDY-EXACT-SIGN)

```
   H_vv = vᵀHess V v = (−D²)_vv + (2ad F_B)_vv + Φ″_vv ,        H_geo := Hess S_W = −D²+2ad F_B .
```

Here `−D²=−D_A^*D_A` is the covariant Laplacian on transverse 1-forms in the adjoint (`D=D(Ā)`), `2ad F_B=[F_B,·]` acts as `ad` on colour and `F` on Lorentz, and `Φ″` is the Haar-log-Jacobian Hessian. The soft `k→0` content is the `k_∥=0`, transverse-magnetic sector.

**Term-by-term sign.**

| term | pointwise sign | expectation | status |
|---|---|---|---|
| `(−D²)_vv=‖D_A v‖²` | `≥0` (Laplacian PSD, adjoint/transverse, gauge-invariant) | `E[(−D²)_vv]≥0` | **真定理**, any measure |
| `(2adF_B)_vv` | **indefinite**; on constant abelian `λ_{n,s}=k_∥²+(2n+1)g|B|+2s\,gB`, bottom `−g|B|` | signed, no sign theorem (below) | pointwise indefinite, **真** (Savvidy) |
| `Φ″_vv` | `≥c_Φ=N/12>0` (SU2 `1/6`, SU3 `1/4`) pointwise | `E[Φ″]≥0`, measured `0.34/0.46` | **真定理** (Lemma B, Haar convexity) |

The only signed term is `2adF_B`. Its colour trace vanishes (`tr_N ad F_B=0`), so `Tr H_{geo}=4Tr(−D²)≥0` (trace-half-positivity), but the operator itself has the Landau-ladder bottom `−g|B|` (depth `∝β|F|`), negative count `N_neg/V=c_N(gB)²/8π²` (`c_SU2=1,c_SU3=3`), **unbounded** as `|F|` grows (PROOF-SAVVIDY-EXACT-SIGN).

### §1.1 The exact Stein / Schwinger–Dyson identity — and the mislabel to correct

Let `f:=∂_nV=⟨v,∇V⟩` be the force along `v`, `H_vv=vᵀHess V v`. On the compact fiber (periodic, no boundary) integrate by parts with `g=f`:
```
   ∫ ∂_n(f e^{−V}) dX = 0  ⇒  E_μ[H_vv] = E_μ[f²]  ⪰ 0 .          (SD / Stein — exact, all β, no convexity)
```
This is the **full** Hessian: `H_vv=(−D²+2adF_B+Φ″)_{vv}`. **It is NOT `E[−D²]=E[f²]`.** `−D²` is only the flat/kinetic summand; `f²=(∂_nV)²` carries the action force *including* `2adF_B` and the `Φ′` cross piece. Equivalently,
```
   E[(2adF_B)_vv] = E[f²] − E[(−D²)_vv] − E[Φ″_vv] .                (budget)
```
Since `E[f²]≥0`, `E[(−D²)_vv]≥0`, `E[Φ″_vv]≥0`:
```
   E[(2adF_B)_vv] ≥ −( E[(−D²)_vv] + E[Φ″_vv] )  =: −budget(v) ,     (budget theorem — 真定理)
   ⟺  E[H_vv] = E[f²] ≥ 0 .                                         (expected Hessian is PSD)
```
So the Savvidy term, **in expectation**, provably cannot drive the expected Hessian below `0`: the budget is never eaten. This is the precise, provable core of "expectation-positive though pointwise indefinite".

> **Correction recorded.** The corpus shorthand `E[−D²]=E[f²]` (and the task's parenthetical) is a **mislabel**: the identity binds the **whole** Hessian `−D²+2adF_B+Φ″` to `f²`, not the Laplacian piece. Read literally it would say the Laplacian expectation alone equals the force energy, forcing `E[2adF_B]+E[Φ″]=0`, which is false (`E[Φ″]≈0.34>0`). The correct object — and the one the numerics/Stein test (`E[H]=3.834`, `E[f²]=3.833`) actually check — is `E[−D²+2adF_B+Φ″]=E[f²]`. This correction matters: it is precisely the deleted `2adF_B` that carries the pointwise negativity and the deleted `Φ″` that is the provable positive budget.

> **Acid test (mandatory).** Free massless `k=0`: `H_vv=A_vv=0`, so `E[H]=E[f²]=0` **and** `Var(a_{k→0})=1/(A⁻¹)_{vv}=∞`. The identity is satisfied with **zero** while the target diverges: it is **gap-blind** (SUB1 §8, B4 §3). It correctly returns `0≥0` on the free field (not spurious positivity) but therefore **cannot be the positivity that fires** — it holds identically on free and interacting theories.

---

## §2 · `E[2adF_B]` at typical field: sign, magnitude, β-dependence

**Pointwise (recap).** `2adF_B` is indefinite with bottom `−g|B|` and `21.7%` negative directions (SU2, W6 `min=−1.12β=−2.81`), `42%` (SU3 β=2.5); after the exact gauge quotient the transverse Hessian is *still* `7.1%` negative at the SU2 crossover (W4, `min −0.82β`) — the negativity is **physical**, not pure gauge (every negative mode is transverse-magnetic with `δF≠0`; PROOF-SAVVIDY-WGEO-DOMINATION §1). It is monotone and unbounded in `|F|` (count `∝∫|F|²`, weight `≲∫|F|³`, depth `−g|B|`); BPS (SD/ASD) backgrounds are PSD, and the equal mixture `Q=0` (Savvidy) is the **most** unstable.

**In expectation (the crux).** The budget identity fixes `E[2adF_B]` exactly as the *residual* of three non-negative quantities:
```
   E[(2adF_B)_vv] = E[f²]_vv  −  E[(−D²)_vv]  −  E[Φ″_vv] .
```
Numerically at the crossover (W7, β-units): `E[H]=E[f²]≈2.45` (SU2) / `3.19` (SU3 β=5.7), `E[Φ″]≈0.34/0.46`, `E[−D²]=O(1)` (the per-mode covariant kinetic). Hence `E[2adF_B]≈2.1−E[−D²]` — of **undetermined sign in principle**, positive if the Laplacian expectation is subdominant. There is **no sign theorem for `E[2adF_B]` alone**, only the **floor** `≥−budget` from the budget theorem. The floor scales like `−E[−D²] ∝ β` at large β (the Savvidy negativity grows with the field), so at strong coupling the margin between `E[f²]` and `budget` can be tight — which is exactly why the unconditioned identity alone does not close anything.

**β dependence.** At `β→0` (pure Haar): `2adF_B→0`, `H→Φ″` only, `E[H]=E[Φ″]≥N/12` — the β=0 closure is exact (`V_eff″=0.248` vs `E[Φ″]−Var(Φ′)=0.246`, PROOF-JACOBIAN-SAVES §2). As β grows, both `E[f²]=E[H]` and `budget` grow; the budget theorem holds at every β but its slack `E[f²]` and the margin `E[f²]−(E[−D²]+E[Φ″])` are not uniformly bounded below by a positive `λ` in any proven way.

**Rigorous lower bound on `E[H]` from `E[f²]`-condensate + `W_geo`?**
- `E[H]=E[f²]≥0` is **rigorous** and needs no condensate input; it is the sharpest *uniform* bound available.
- A stronger `E[H]≥λ` would require a proven **lower bound `E[f²]≥λ`** on the force/condensate second moment. At `β=2.5` this is unproven; the closest proven positivity is the `β<1/48` residue `f₁²>0` (P3), which is the **atom residue**, not `E[f²]`.
- **Do not add `W_geo` to `E[H]`.** `W_geo=E[Φ″]−Var(Φ′)≥0` is the curvature of the *conditioned* `Φ`-marginal (a Prékopa average), whereas `E[H]=E[f²]` is the *unconditioned* full-Hessian identity; `Φ″` is already inside the Stein identity. Summing `E[f²]` with `W_geo` is a category mix across conditioning levels. The rigorous, level-correct bound is simply `E[H]≥0`.

---

## §3 · The sharpest provable conditional — and why it is a variance bound, not a sign

**The task's proposed conditional is not a reduction.** "If `E[2adF_B]≥−(budget−λ)` then `G″≥λ`": by the budget theorem the antecedent is **exactly** `E[f²]≥λ` (the *unconditioned* force energy clears `λ`). But `G″(0)=E[H|0]−Var(f|0)+W` is a **conditioned** quantity. The unconditioned `E[H]≥λ` does **not** imply it. A clean counterexample: on any measure, `G″(0)=−(log p₀)″(0)` where `p₀` is the marginal of `a_{k→0}`. Take `p₀` even and **bimodal** (a dip at `0`): `p₀″(0)>0⇒G″(0)<0`, while `E[H]=E[f²]≥0` still holds identically. So `E[H]≥0` (even `≥λ`) is consistent with **any** sign of `G″(0)`. The proposed implication is **false**: the unconditioned expected-Hessian sign does not control the conditioned marginal curvature.

**The correct sharpest conditionals** (each an *exact* identity, PROVE-Gpp §7, PROOF-JACOBIAN-SAVES §1, INVENT-VARIANCE-INEQ §1):
```
   G″(0) = E[H_vv|a=0] − Var(f|a=0) + W_full        (W_full = W_geo+X ≥ N/12)            [Jacobian split]
         = S(v) − D(v)                              (S = Schur/1/Var(a), D = fiber BL defect)  [SCBL]
         = −(log p₀)″(0) ,     Ward: G″·Var(a)=1 ⟺ marginal Gaussian                    [marginal]
```
Sufficient forms:
- **(C1) SCBL:** `D(v) < S(v) − λ  ⇒  G″≥λ`. Measured `D/S≈2%` (`S≈1/Var(a)≈0.27–0.45`), so `G″≈0.98 S>0` — but `D<S` uniform is the fiber BL defect `D=Var(f|0)−E[⟨K⁻¹c,c⟩|0]`, finite iff the indefinite fiber `K=Π₀HessVΠ₀` has a spectral gap. **This *is* H1.**
- **(C2) Jacobian:** `Var(f|0) − E[H_vv|0] < W_full − λ  ⇒  G″≥λ`. Measured `Var−E[H]∈[−0.06,+0.02]≈0 < W_full≈0.22–0.31`, margin `≈0.3`. Proving the uniform upper bound on `Var(f|0)` — **this *is* H1.**
- **(C3) Ward:** `G″≥λ ⇔ Var_μ(a_{k→0})≤1/λ` (one direction is a theorem via Brascamp–Lieb; the converse holds iff the marginal is Gaussian). The zero-momentum susceptibility bound — **this *is* H1.**

All three are one object: a **one-sided upper bound on the conditioned force variance** `Var(f|0)` (= `D`/`Var(a)`). **The uncontrolled piece is a VARIANCE, not a sign.** The expectation-level sign of `E[H]` is *already settled* (provably `≥0`, `≈+2.45`); what remains is the *conditioned variance*.

**Does `D/S≈2%` already evidence expectation-positive curvature?** **Yes, with a precise reading.** `D/S≈2%` says the coarse marginal is `98%` Gaussian, i.e. `G″=S−D≈0.98/Var(a)>0` — the **expectation-level coarse curvature is positive with a `2%` margin**, in agreement with `E[H]≈Var≈2.45` (`E[H]−Var≈0`) and `G″_{meas}=+0.27,+0.32,+0.51`. So the *expected* Hessian's positivity is not merely `≥0`; it is `>0` with margin at the measured β. **But** `D<S` is *measured*, and `D` depends on the inverse fiber Hessian `K⁻¹`, which blows up precisely at the edge where a transverse massless mode crosses (`det K=0 ⇒ ‖K⁻¹‖→∞`). Proving `D<S` uniformly is the fiber spectral gap = H1, **not** a sign statement about `E[2adF_B]`.

---

## §4 · Pointwise indefinite vs expectation-positive — precisely where the gap lives

Four levels, in increasing conditioning:

| level | object | sign | status |
|---|---|---|---|
| **L1 microscopic, pointwise** | `H_vv` on all of `TM` | **indefinite**, bottom `−g|B|=−1.12β`, `21.7%` neg (W6), `7.1%` after gauge quotient (W4) | **真** (Savvidy, physical) |
| **L2 coarse, pointwise** | `PᵀH_vvP` (block/collective) | positive at crossover (`+0.4β` W6, `+0.81β` blocked W4, `0%` neg), SU3 deep-SC still `12%` neg | pointwise-over-samples, not uniform |
| **L3 expectation, unconditioned** | `E[H_vv]=E[f²]` | **`⪰0` PSD**, `≈+2.45` | **真定理** (Stein/SD), but **gap-blind** |
| **L4 gap, conditioned** | `G″(0)=E[H|0]−Var(f|0)+W` | measured `+0.27…+0.51`, carried by `W` | **=H1** |

**The gap lives between L3 and L4, not between L1 and L3.** At the `P`-symmetric vacuum `E f=0`, so unconditionally `Var(f)=E[f²]=E[H]` **exactly** (equipartition) — the flat object `E[H]−Var(f)` is **identically `0`**, and the entire `G″` is the **conditioning correction** `W_geo+X` plus the change in slice variance. The pointwise Savvidy negativity (L1) is **real**, but its expectation is *exactly compensated* by the force second moment `E[f²]` (Stein, L3); the compensation is **unconditioned** and does **not** transfer to the `a=0` slice. "Pointwise indefinite ⇒ expectation-positive" is thus a **true and proven** transition — but it crosses only from L1 to L3, and L3 is **on the free-field side of the gap** (`E[H]=E[f²]≥0` holds for the free massless field too, with `E[H]=0, Var(a)=∞`). The gap is the further **conditioning** L3→L4 (from the equilibrium measure to a fixed coarse slice), whose sign is the non-Gaussianity of the coarse marginal at typical action (ABC-NEW-INEQ §1; Ward `G″·Var(a)=1` iff Gaussian).

So the honest slogan: **the pointwise Savvidy wall is compensated in expectation by a true identity; the gap is the one step that the identity cannot see — conditioning.** `E[2adF_B]` neither eats the budget (budget theorem) nor signs the gap (conditioning, not sign).

---

## §5 · Scorecard (5-point, corpus convention)

1. **Nonperturbative — PASS.** Budget theorem / Stein identity exact at every β, no `g`/`1/N`.
2. **Sign-controlling — PARTIAL (expected-Hessian sign only).** Proves `E[H]⪰0` and hence `E[2adF_B]≥−budget`; does **not** prove `G″≥λ` (the conditioned sign), which is a variance bound.
3. **Non-circular — FAIL as a gap tool.** Closing `G″≥λ` from here forces `Var(f|0)≤E[H|0]+W−λ` (=H1); the unconditioned identity is silent on `Var(a)`.
4. **Non-SUSY + fails free field — FAIL (vacuous / gap-blind).** Free massless saturates `E[H]=E[f²]=0` while `Var(a)=∞`; the bound is met by the free field, so it does not fail the acid test in the separating sense.
5. **Uniform `m_eff²(k→0)≥λ` — FAIL.** Only `E[H]≥0` uniform; `E[H]≥λ` needs the unproven condensate bound `E[f²]≥λ`, and `G″≥λ` needs H1.

**Net: 1 PASS, 1 PARTIAL, 3 FAIL — all three failures are the single conditioning/variance species.** A genuine theorem about the expected Hessian (`L3`), not the gap.

---

## §6 · Sharpest line

The exact Stein identity `E_μ[Hess V]=E_μ[∇V⊗∇V]` makes the **full expected Hessian provably PSD** at the soft `k→0` mode — `E[−D²+2adF_B+Φ″]=E[f²]⪰0`, hence `E[2adF_B]≥−(E[−D²]+E[Φ″])` **exactly** — so the Savvidy term genuinely does **not** eat the budget and "expectation-positive though pointwise indefinite" is **true and closed**. But it is the **unconditioned** identity, and at the `P`-symmetric vacuum it reads `E[H]=Var(f)` **identically** (the flat object is exactly `0`), so the entire gap is the **conditioning correction** `D<S` / `Var(f|0)≤E[H|0]+W−λ` — a one-sided **variance** bound (=H1), not a sign of `E[2adF_B]`. The gap sits **between L3 (provable expectation-positivity) and L4 (conditioned curvature)**, both sides of which are on the wrong side of the free-field acid test (`L3` is saturated by the free field). The measured `D/S≈2%` evidences `G″≈0.98 S>0` — expectation-positive with margin — but its proof is the fiber spectral gap, i.e. H1 again.

---

## §7 · The exact stuck step

Prove, uniformly in `(L,a)` and for the physical coarse mode `v⊥orbit`:
```
   (H1-cond)   Var_μ(∂_nV | a=0) − E_μ[H_vv | a=0]  <  W_full(v) − λ ,
```
equivalently `D(v) < S(v) − λ` (SCBL) or `Var_μ(a_{k→0}) ≤ 1/λ` (Ward). **One object, three faces: the conditioned force-variance ceiling.** Everything on the *expectation* side is already closed (`E[H]=E[f²]⪰0`, budget theorem); everything on the *pointwise* side is correctly indefinite (Savvidy); the single open step is the **conditioning** that bridges them — and it is a variance bound, not a sign. No sign of `E[2adF_B]` — proven, measured, or conjectured — supplies it.

*Refs.* Savvidy, Phys. Lett. B **71** (1977) 133; Nielsen–Olesen, Nucl. Phys. B **144** (1978) 376; Leutwyler, Nucl. Phys. B **179** (1981) 129; Bogomol'nyi, Sov. J. Nucl. Phys. **24** (1976) 449; Belavin–Polyakov–Schwarz–Tyupkin, Phys. Lett. B **59** (1975) 85; 't Hooft, Phys. Rev. D **14** (1976) 3432; Brascamp–Lieb, Adv. Math. **22** (1976) 366; Prékopa, Acta Sci. Math. **32** (1971) 301; Efron–Stein (1981); Stein's lemma / Schwinger–Dyson IBP. Internal: PROOF-SAVVIDY-EXACT-SIGN (`−D²+2adF`, ladder `−g|B|`, count/depth); PROOF-SELFDUAL-BACKGROUND (BPS PSD, `Q=0` most unstable); PROOF-SAVVIDY-WGEO-DOMINATION (negative fiber physical after quotient); PROOF-JACOBIAN-SAVES (`G″=Δ_dyn+W_geo+X`, `W_geo≥N/12`, β=0 closure 0.248/0.246); FORGE5-CONE-VARIANCE (§"exact identity" `E[H]=E[f²]`, gap-blind); SUB1-COV-BOUND (§3–8, equipartition `Var(f)=E[H]`, free-field `E[H]=E[f²]=0, Var(a)=∞`); INVENT-VARIANCE-INEQ (SCBL `G″=S−D`, `D/S≈2%`); PROVE-Gpp-positive (§7 exact `G″=E[H|0]−Var(G|0)+W`, §3 `−Cg|F|` floor); B4-VARIANCE-MAG (`VarG≈2` L-stable, `E[H]−Var∈[−0.06,+0.02]`); W6-BE-at-2.5 (`min −1.12β`, 21.7%/42% neg); W4-blocked-positivity (`7.1%` neg after `Q⊥`, `+0.81β` blocked); W7-Gpp0 (`G″` +0.27/+0.32/+0.51, `E[H]≈Var≈2.45`, `1/Var≈0.27`); PROOF-COARSE-POINCARE (H1, `PᵀHP` min −0.3); OTHER-SAVVIDY-SIGN (parity `Q=0` balance); REDUCTION-AUDIT (H1 = susceptibility, one-way).
