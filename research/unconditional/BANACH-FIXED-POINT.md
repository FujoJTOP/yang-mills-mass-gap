# BANACH-FIXED-POINT — constrained curvature as a self-consistent Riccati fixed point: Schauder banks β=0; the zero mode is neutral at k→0; Savvidy makes T expansive on the shell R_k must absorb

> **2026-11-23 · NEW SPECIES (nonlinear analysis / Banach fixed point).** 判词. Write `κ(a)=V_eff″(a)` — via `r=−½V_eff′` — as the fixed point of a self-consistent Riccati/constrained-partition map `T` on a weighted Hölder space, and ask whether contraction/Schauder/monotone analysis forces `κ≥λ` uniformly. **The vehicle is genuine; β=0 becomes a real fixed-point theorem. The target fails at 2.5 structurally: `T` is a parabolic (high-frequency) contraction whose zero-frequency direction is neutral at `k→0`, and the Savvidy negative mode makes `T` expansive exactly where `R_k` must absorb it.** The acid test passes — but as the trivial degenerate fixed point a contraction near `κ=0` *certifies* (PROOF-FRG-DSE-GAP-EQUATION: "proves the trivial fixed point is unique — the wrong answer"). **#1,#4 PASS; #2(2.5),#3,#5 FAIL. No closure.**

## §1 The equation
`V_eff=−2 log u`, `u=p^{1/2}`; `r=−½V_eff′=u′/u` obeys the Riccati equation `r′+r²=W−E`, `W(a)` the Bakry–Émery potential of the `a`-fiber, and `κ=V_eff″=−2r′`. **Uniform positive curvature ⟺ uniform Riccati decay `r′≤−λ/2` ⟺ the zero-mode drift `b(a)=G″a` contracts (`G″≥λ`).** Self-consistency: `W` is taken over the constrained measure, hence depends on `r`, giving the nonlinear map `T: r↦soln(r′+r²=W[r]−E)` on `B=C^{2,α}_w([−A,A])` (Gaussian weight), `T` compact on bounded sets.

## §2 The contraction constant
`DT` is the linearized constrained second variation — the regularized propagator `(Γ_k^{(2)}+R_k)⁻¹(k∂_kR_k)(Γ_k^{(2)}+R_k)⁻¹` (E1-FRG §1). On every shell `k>0`, `Γ_k^{(2)}+R_k⪰0` (Dupuis et al.), so `ρ(DT)<1`: contraction holds. But it is *frequency-localized*: the `q→0` zero mode, which *is* `κ(0)`, is neutral, and as `k→0`, `R_k→0` gives `ρ(DT)→1`. The machine hands back exactly the residual curvature it set out to bound (E1's vanishing-regulator verdict): stable for `k>0`, neutral at `k=0`.

## §3 Does Savvidy break it? Yes
The Savvidy Hessian `−D²+2adF_B` has transverse-magnetic eigenvalue `−2g|F_B|` on `|k|<√(2g|F_B|)`. There `(Γ^{(2)}+R)` is dominated by `R_k` only until `R_k≈2g|F_B|`; the factor `(k∂_kR)/(R−2g|F_B|)²` **blows up** as `R_k↓2g|F_B|`, so `sup_fiber ρ(DT)=∞` — **no uniform contraction, and `{κ≥λ}` is not `T`-invariant**: a `κ≥λ` profile maps to `κ<0` on the typical-field tail (D3: log-concavity fails at 2.5). Schauder/Leray–Schauder give *existence* in a ball but **no sign, no uniqueness** — the ball giving `κ≥λ` is the non-invariant cone. Monotone-operator (Krasnosel'skii) needs an invariant order cone: SU(N) has none (torsion, C8 §b), and `Cov(G,Φ′)>0` (C8) is order-*reversing*. Riccati comparison against `W₀≥N/12` fails because `W` is indefinite on the cone — the ODE face of PROOF-OPTIMAL-CONTROL-HJB's "indefinite `Hess S_W` has no stabilizing solution."

## §4 Where it is a theorem — β=0
At β=0, `W=Φ″⪰N/12` is convex; Prékopa makes the marginal log-concave, so `{κ≥N/12}` *is* invariant and `T` order-preserving: Schauder/Krasnosel'skii give `κ*≥N/12`, uniform in `L` — the β=0 gap as a **fixed point**, not a point value (upgrade of k0-CONSTRAINT-EFFPOT).

## §5 Acid test
Free massless: `W=0`, `r≡0`, `κ=0`, `T(0)=0`; the contraction near `0` certifies the trivial fixed point. A nonzero branch needs an eigenvalue of `DT` crossing `1` (bifurcation) — the Schwinger pole — whose insertion *is* `f₁²>0` (circular). Acid test PASSED, vacuously: the only contractively-certified fixed point is the gapless one.

## Five-point scorecard
1. **Nonperturbative — PASS.** Exact Riccati + Schauder; no expansion.
2. **Sign lower bound `κ(0)≥λ` — FAIL at 2.5.** Cone `{κ≥λ}` not `T`-invariant; β=0 gives `N/12`.
3. **Non-circular — PASS as vehicle, FAIL at closing.** The 2.5 fixed point needs `Var(G|0)` = H1; the nonzero branch is the bifurcation = the order parameter.
4. **Non-SUSY + fails free field — PASS (clean).** `W=0⇒κ=0`; Savvidy cone non-SUSY.
5. **Uniform `≥λ` in (L,a) — FAIL.** Only β=0 `N/12`; zero mode neutral at `k→0`.

**Sharpest line.** Nonlinear analysis is the right *calculus* — the constrained curvature is a self-consistent Riccati variable, β=0 is a genuine Schauder–Krasnosel'skii fixed point with invariant cone `κ≥N/12`, and the free field collapses to `κ=0`. But the operator is a high-frequency contraction whose zero mode is neutral at `k→0`, and Savvidy makes it expansive on the exact shell `R_k` must absorb: "uniform positive curvature" lives on the **neutral** direction of `T`, unreachable by contraction. **Savvidy breaks the contraction; the machinery certifies the wrong (gapless) fixed point or returns the `R_k→0` curvature = H1. #1,#4 pass; #2(2.5),#3,#5 fail — a sharper name for the wall, not a crossing.**

*Refs.* Banach, *Fund. Math.* **3** (1922); Schauder, *Studia Math.* **2** (1930); Leray–Schauder, *Ann. Sci. ENS* **51** (1934); Krasnosel'skii, *Positive Solutions of Operator Equations* (1964); Reid, *Riccati Differential Equations* (1972); Bakry–Émery, *LNM* **1123** (1985); Savvidy, *Phys. Lett. B* **71** (1977). Internal: E1-FRG (convexity, vanishing regulator); k0-CONSTRAINT-EFFPOT (`V_eff″=G″`, β=0 `N/12`); D3 (log-concavity fails at 2.5); C8 (torsion, `Cov>0`); F1 (barrier = H1); F2 (rate = H1); PROOF-FRG-DSE-GAP-EQUATION (trivial fixed point); PROOF-OPTIMAL-CONTROL-HJB (no stabilizing Riccati); PROOF-COARSE-POINCARE (H1); PROOF-JACOBIAN-SAVES.
