# C-BBS-SUFFICIENCY — the capacity/conductance floor is the Poincaré constant of the *coarse marginal* = `1/χ = G″` = **H1** (necessary), not the physical spectral gap `Δ = inf supp ν` (sufficient); C-BBS cannot deliver (MS)

> **2026-11-23 · crux ruling.** Target: does the C-BBS capacity/Dirichlet trap (UV trapped `g_j→0`; thermal-shell Arrhenius `O(1)`) give the **sufficient** `(MS) ⇔ inf supp ν>0`, or again only the **necessary** `H1 ⇔ Var(a)≤1/λ`? **判词 (verdict): NECESSARY, H1.** The conductance floor `Φ_eq≥c` yields `λ_1 ≥ Φ_eq²/4` (Sinclair), but `λ_1` is the **Poincaré constant of the 1D coarse marginal** `ν(da)=e^{−V̄}da/Z`. Tested on the single coordinate `a` (unit carré du champ `Γ(a)≡1`), Poincaré reads `Var_ν(a) ≤ 1/λ_1`, equivalently `λ_1 ≤ 1/Var(a)=1/χ`. A **lower** bound `λ_1≥c` therefore gives `χ ≤ 1/c`: a bounded **susceptibility**, verbatim H1. `λ_1` is the **spectral gap of the effective marginal diffusion**, which equals `G″=1/χ` (Gaussian-exact, `G″χ≈1`), *not* the transfer-matrix gap `Δ=inf spec(H_phys)^{1/2}`. Sufficiency would need the Poincaré constant of the **full gauge-quotiented measure** (all test functions `f`), and that is exactly where the Savvidy indefinite Hessian returns and the capacity dictionary dies. **C-BBS can at best close a uniform H1; the edge needs the second-level isolation `f₁²>0 = H*`, invisible to every bottom/Dirichlet functional.**

## §1 · The dictionary, followed to the end
Capacity replaces "sign an indefinite Hessian" by "minimize a convex Dirichlet energy," and the sign is free (`λ_1≥Φ_eq²/4⩾0`). But the object it bounds is the **marginal** `ν`. The two exact facts:
```
Var_ν(a) ≤ (1/λ_1)·ℰ(a,a),   ℰ(a,a)=∫Γ(a)dν=1   ⇒   λ_1 ≤ 1/χ = 1/Var(a)      [Rayleigh, one f]
λ_1 := inf_{Var>0} ℰ(f,f)/Var_ν(f)                                        [Poincaré constant of ν]
capacity floor  λ_1 ≥ Φ_eq²/4 = c  ⇒  χ ≤ 4/Φ_eq².                        [= H1, upper bound on χ]
```
So the floor is **an upper bound on the variance `χ`** — precisely the necessary condition, nothing more. POTENTIAL-THEORY-CAPACITY itself scores #5 FAIL ("only a bootstrap on `λ` … the `L`-uniform floor `Φ_eq≥c` is exactly H1"); this note makes the *reason* rigorous and answers the spectral-vs-susceptibility crux.

## §2 · The crux: Poincaré constant ≠ spectral edge
A Poincaré constant of the **full physical measure** `μ` (bounding *all* gauge-invariant `f`) *is* sufficient: `Var_μ(f)≤λ_1^{-1}ℰ(f,f)` ⇔ exponential decay of connected correlators ⇔ gap (Shen–Zhu–Zhu; G2H (P1)). The capacity construction does **not** bound that object. It tests **one** function, the zero-momentum coordinate `a`, against the **1D marginal** `ν`, and a single-function Rayleigh quotient only pins `λ_1` from *above* (`λ_1≤1/χ`). Two structural reasons the edge is lost:
- **Wrong measure.** `ν` is the *equal-time* law of one mode. By Feller, *every* probability on ℝ is the invariant measure of some gapped 1D diffusion, so `ν`'s best Poincaré constant is a functional of a single-time distribution — it cannot encode the spacetime spectral edge `inf supp ν_O(s)`.
- **Wrong functional class.** The physical mass gap is `Δ=inf supp ν_O`, a **support** datum. Poincaré is an **L²/second-moment** datum. They differ by the residue/atom weight: `χ=∫s^{-1}dν_O(s)` is a residue-weighted *harmonic* mean of masses, dominated by the lightest **residue**, not the lightest **mass**. In the Gaussian (`J=S=G″`) limit `λ_1=G″=1/χ` literally **is** the susceptibility; off-Gaussian, Jensen makes the harmonic-mean stiffness *strictly stronger* than the centre `G″(0)`, yet still an H1-type second-moment bound (D3 §4).

## §3 · The susceptibility trap, reinstated
The audit counterexample survives the capacity translation. Take `ν_O=wδ_ε+(1−w)δ_S`, `w/ε≤1/λ`. Then `χ=w/ε+(1−w)/S<∞`, so the marginal has a **positive** Poincaré/conductance gap (`λ_1>0`, `Cap>0`, finite `Φ_eq`), while `inf supp ν_O=ε→0`: **no gap.** A capacity floor, uniform or not, cannot rule this out — it sees `χ`, not the thin `δ_ε` mass (G2H §2: moments/susceptibility do not imply support). This is the exact mirror of E2's finding that information/capacity ceilings bound `Var(a)` **from below** (rate–distortion) or only the conditional variance at `snr>0` — never the prior from above without already containing `χ`.

## §4 · The constructive C-BBS idea: UV trapping ⊥ the IR edge
"UV trapped (`g_j→0`)" is the BBS Gaussian-anchor side — orthogonal to the IR edge (P2-BBS-to-gauge: the untransplantable hypothesis is the Savvidy-safe large-field bound). "Thermal-shell Arrhenius `O(1)`" removes the exponential barrier at the target scale, so a bootstrap on `λ` can in principle be closed — **delivering a uniform lower bound on `λ_1`, i.e. a uniform H1.** But uniform H1 is still necessary-not-sufficient (G2H (P2); REDUCTION-AUDIT (b)). Closing (MS) additionally needs `f₁²>0` with isolation `I2` (`H*`), a **two-level** (`λ₁−λ₂`) datum that every bottom/Dirichlet/Poincaré/Bushell functional is blind to (PROOF-EDGE-E2; G2G §4). The capacity floor is one rung too low.

## Scorecard
1. **Nonperturbative — PASS.** Capacity/Wiener/max-flow exact, all β.
2. **Sign-controlling `G″≥λ` — PASS(shape, `≥0`), FAIL(uniform).** Convex/positive min, but only certifies `χ`, not `Δ`.
3. **Non-circular — FAIL.** Uniform `Φ_eq≥c` = uniform H1 in flow costume.
4. **Non-SUSY + fails free field — PASS (cleanest).** Flat free zero mode ⇒ `Cap=0`.
5. **Uniform `inf supp>0` (MS) — FAIL.** Delivers bounded `χ`; the edge/isolation is not a bottom functional.

## Exact failing estimate
```
λ_1 ≥ Φ_eq²/4  [capacity floor]      λ_1 ≤ 1/Var(a)=1/χ   [Rayleigh on a]   ⇒  χ ≤ 4/Φ_eq²   [H1]
Gaussian:  λ_1 = G″ = 1/χ  ⇒  capacity floor ⟺ G″(k→0)≥c ⟺ H1            [same number]
(MS) needs: inf supp ν_O ≥ Δ²  ← support edge, NOT second moment;  §3 counterexample defeats any χ-bound.
free massless: V̄≡0 ⇒ Cap=0 ⇒ λ_1=0 ⇒ χ=∞  [acid PASS, but the content is χ]
sufficiency needs Poincaré of FULL μ (all f) ⇒ Savvidy indefinite ⇒ capacity dictionary dies.
missing:  f₁²>0 + isolation I2  = H*  = two-mass λ₁−λ₂  [EDGE-E2 bottom-blind]
```

**Sharpest line.** The capacity floor is the **cleanest non-vacuous acid test** in the species family (flat free zero mode ⇒ zero capacity) and it buys a genuinely convex, sign-free route to a **bounded susceptibility** — but a Poincaré constant is the spectral gap of the *marginal diffusion*, and for that measure it **is** `1/χ=G″`, the H1 number; it is not the physical spectral edge `Δ=inf supp ν`. The "spectral gap ⇒ mass gap" step is the fallacy: it holds for the full gauge-quotiented measure and fails for the 1D marginal, where sufficiency is lost in the reduction to a single test function. C-BBS (`g_j→0` UV + `O(1)` thermal shell) can at best prove a **uniform H1**; `(MS)` still needs the second-level isolation `f₁²>0`, which no bottom capacity functional can see. **Same lock, new (and sharper) key shape. Necessary, not sufficient.**

*Refs.* POTENTIAL-THEORY-CAPACITY (verdict: `Φ_eq≥c` = H1; #5 FAIL); REDUCTION-AUDIT (chain: susceptibility necessary, `H*` sufficient; `ν=wδ_ε+(1−w)δ_S` counterexample); G2H-BLIND-CLEAN §2–5 ((P2) moments≠support; (MS)); G2G-TRANSFER-SPECTRAL §4 (`λ₁−λ₂` bottom-blind; `G″=H1`); D3-SUSCEPTIBILITY-VARIATIONAL (`1/Var` vs `G″`; Jensen); E2-INFORMATION §2–4 (capacity/entropy bounds `Var` from below / at `snr>0`); CONSTRUCTION-METHODS (a)+(c) BBS trap; P2-BBS-to-gauge (Savvidy-safe LF untransplantable); PROOF-EDGE-E2 (Temple/Bushell regress). External: Sinclair, *Combin. Probab. Comput.* **1** (1992) 351; Jerrum–Sinclair, *SIAM J. Comput.* **18** (1989) 1149; Feller, *Ann. of Math.* **55** (1952) 468; Wiener, *J. Math. Phys.* **3** (1924) 127; Kakutani (1944); Doyle–Snell, Carus **22** (1984); Brascamp–Lieb (1976); Bauerschmidt–Brydges–Slade, *LNM* **2242** (2019).
