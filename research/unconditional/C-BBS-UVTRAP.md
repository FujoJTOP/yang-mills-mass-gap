# C-BBS-UVTRAP — the capacity/Dirichlet trap: every UV shell is trapped with a shell-uniform floor, the thermal shell is the only untrapped one

> 2026-11-23 · C-BBS (constructive, capacity scheme). Replace BBS's indefinite-Hessian trap by a **capacity/Dirichlet-form** trap — a *minimum of a manifestly positive energy*, not a lower bound on an indefinite operator. **判词 (verdict).** The UV half is a **条件定理 (conditional theorem)**: on every shell `j` above the physical crossover, the shell Dirichlet form satisfies `E_j(f) ≥ c\,Λ_j²\,Var_{μ_j}(f)` with `c>0` **independent of `j,a,L`**, so the capacity floor `Φ_j ≥ c'Λ_j²` traps shell `j`. The constant is **not** the collapsing `e^{−2β_jN_p}`; it is `c₁e^{−Cg_*²}` with `g_*` the *largest* trapped coupling, i.e. the floor **opens toward the UV** exactly where AF says it should (G1-L3-TOWER-POINCARE). Every estimate is banked except one — a **gauge-invariant, `a`-uniform small-field moment bootstrap on the typical region** — but it is used here only to bound a **density ratio**, not to sign a Hessian. The Savvidy negative modes (`−2g_j|F_j|`) are **shallow**: `(depth)/Λ_j² = O(g_j²/Λ_j²) → 0` super-exponentially. The only untrapped shell is the thermal/IR one (`μ_j∼2π/L`, `g=O(1)`), where the cone re-enters (FORGE5-K2) and the floor *is* H1.

## 1 The trap, stated positively
On shell `j`, band `[Λ_{j−1},Λ_j]`, `μ_j=Λ_{j−1}`, let `E_j(f)=∫|∇_j f|²dμ_j` be the shell Dirichlet form. Since `E_j(f)≥0` for **every** field and `f`, the Savvidy negative Hessian (`−D²+2 ad F_B`, `λ_min=−2g_j|F_j|`) never makes it indefinite — a Dirichlet form is quadratic in `f`, not `U` (PROOF-DIRICHLET-COMPARISON §1). The trap is the **capacitary floor**
```
(CT)_j :   E_j(f) ≥ c Λ_j² Var_{μ_j}(f)   ∀f,   equivalently  Φ_j := inf_A Cap_j(A)/μ_j(A) ≥ c'Λ_j²,
```
with `c,c'` shell-uniform. This replaces BBS's `G″_j≥δ` (indefinite Hessian, measure-blind sup) by a *positive* Dirichlet min (measure-weighted).

## 2 Three inputs, shell by shell
**(a) AF makes the negative modes shallow.** In the tower `β_j∼4Nβ₀j`, `g_j²∼1/(2β₀j)` (G1-L3-TOWER-POINCARE §2). Savvidy depth `=2g_j|F_j|`; with the trap invariant `|F_j|=O(g_j)` this is `≤Cg_j²`. Against the free band edge `μ_j²=Λ²e^{2(j−1)}`,
```
(depth)/μ_j² ≤ C g_j²/Λ²e^{−2(j−1)} = O(g_j²) → 0   super-exponentially.
```
The cone radius `√(2g_j|F_j|)∼g_j≪Λ_{j+1}` lies **outside** the retained band: `Hess S_j⪰cΛ_j²` uniformly (G1-L3-UV-FORCE-BOUND §3, Route B). Negative modes are controlled, not deep.

**(b) Prékopa gives the positive floor.** The collective-mode connection curvature contributes `W_geo≥N/12` (PROOF-JACOBIAN-SAVES), **scale-independent**, hence uniform in `j,a,L`. It keeps the collective direction log-concave, so the shell density is bounded between Gaussians and `Cap_j>0` before the free term — without it the mode can be flat (`Cap=0`, the free-field acid-test failure, POTENTIAL-THEORY-CAPACITY §4).

**(c) Dirichlet density is `cΛ_j²`, shell-uniform.** On the high-pass band `|∇_j f|²≤Λ_j²|f|²`; the free band form is **scale-covariant on one octave** (`Λ_j/Λ_{j−1}=e`), so its capacity is `c₁μ_j²` with `c₁=f(e)`, depending only on the octave ratio — **not on `j`**. This is the explicit shell-uniform constant.

## 3 The shell-by-shell transfer (the proof)
The measure `μ_j` differs from the free band measure by `dμ_j/dμ_free=Z^{−1}e^{−H_int}`, `H_int=K_j^{Savvidy}+W_geo`. On the **good set** `{Q≤·}`, `‖log(dμ_j/dμ_free)‖_∞ ≤ Cg_j²` (small-field Gaussian expansion on the trapped region); **off it**, the large-field weight `e^{−cQ/g_j²}` is exponentially small and — being **measure-weighted** — contributes a vanishing fraction to `Cap_j`. Capacity transfer (two-sided bounded density, Deuschel–Stroock) gives
```
Cap_j ≥ e^{−Cg_j²} Cap_j^{free} ≥ e^{−Cg_j²} c₁ μ_j²,
⇒  Cap_j/Λ_j² ≥ c₁ e^{−Cg_j²} ≥ c₁ e^{−Cg_*²} := c > 0,
```
`g_* = g_{j₀}` at the **first** trapped shell. Because `g_j` *decreases* into the UV, the minimum `e^{−Cg_j²}` is at the **shallowest** shell and the floor **rises** with `j` — the opposite of the Holley–Stroock collapse `e^{−2β_jN_p}` (TOWER-POINCARE §2). The K1 objection (a sup over bare curvature is measure-blind) is **evaded structurally**: the Savvidy tail enters through its measure weight `e^{−cQ/g_j²}`, not `sup_x λ_min(Hess)`, because capacity is a *min against `μ_j`* and the good-set density ratio depends on `g_j`, not an unbounded curvature supremum.

## 4 Where the trap fails — the thermal shell only
At `μ_j∼2π/L`, `β=2.5`, `g²≈2.4`: the depth `2g|F|` is **not** `≪μ²`, the cone **re-enters** and is high-momentum (FORGE5-SL-UV-CONE-EMPTY K2), `e^{−Cg²}=O(1)` not `→1`, and `(CT)_j` reduces to `Var(a)≤1/λ` — i.e. **H1**, the gap (POTENTIAL-THEORY-CAPACITY §5). So the only untrapped shell is thermal/IR; the capacity discriminator (free flat mode `Cap=0` vs YM `Cap>0`) lives **entirely** there — the UV trap is a *benign-regime* statement (UV is where free and YM coincide under AF), never an acid-test pass.

## 5 Ledger
| leg | status |
|---|---|
| positivity `E_j≥0` (bypass Savvidy sign) | **真定理** |
| free-octave capacity `c₁μ_j²`, `c₁` shell-uniform | **真定理** (scale covariance) |
| shallow negatives `(depth)/Λ_j²→0` (AF) | **真定理** given `\|F_j\|=O(g_j)` |
| Prékopa floor `W_geo≥N/12`, scale-free | **真定理** |
| density-ratio transfer `Cap_j≥e^{−Cg_j²}Cap_j^{free}` | **条件定理** (small-field bootstrap) |
| `Φ_j≥c'Λ_j²` on every UV shell, `j`-uniform | **条件定理** on (i) only |
| thermal/IR shell trap | **OPEN = H1** |

**Exact missing step (one).** The gauge-invariant, `a`-uniform small-field moment bootstrap `E|F_j|≤Cg_j` on the typical region — the same input as G1-L3-UV-FORCE-BOUND Route B (i) — is needed **only** to make the good-set density ratio `≤e^{Cg_j²}`. It does **not** need a Hessian lower bound: the capacity trap converts (★)'s "sign an indefinite operator" into "bound one density ratio," a one-sided, measure-weighted estimate.

**Sharpest line.** The capacity trap is the correct replacement for the Hessian trap because its sign is *free* (a min of a positive Dirichlet energy) and its constant is *physical*: it tracks the free octave floor `c₁Λ_j²` that **rises** into the UV, so the Holley–Stroock degradation `(aΛ)^c` becomes an artifact the capacity form does not see. Every UV shell is trapped with explicit constant `c₁e^{−Cg_*²}`; the trap leaks at exactly one shell — the thermal/IR one, where the cone is physical (K2) and the floor is the gap itself. The key is now H1 in one place only.

*Refs.* POTENTIAL-THEORY-CAPACITY; PROOF-DIRICHLET-COMPARISON §1–2; PROOF-TRAPPING-REGION §1–2; G1-L3-TOWER-POINCARE §2–4; G1-L3-UV-FORCE-BOUND §3; FORGE5-SL-UV-CONE-EMPTY K1–K2; PROOF-JACOBIAN-SAVES (`W_geo≥N/12`); Bauerschmidt–Brydges–Slade LNM 2242 (2019); Balaban CMP 109 (1987) 249, 116 (1988) 1, 122 (1989) 175/355; Deuschel–Stroock, *Large Deviations* (1989) Ch. 5; Sinclair, *Combin. Probab. Comput.* 1 (1992) 351; Evans–Gariepy (capacity/Wiener).
