# C-BBS-FORMULATE — capacity/Dirichlet-form trap: the indefinite Hessian is never signed; the residual is a per-scale conductance floor (= H1 in flow costume)

> 2026-11-23 · SHAPE-C formula. Successor to FORGE5-CONSTRUCTIVE / CONSTRUCTION-METHODS (c): replace the Hessian-based BBS trapping region — whose sign step `δ²K_j⪰δ` dies on the Savvidy cone (`δ²R_j=O(g²)≪λ_min=O(μ²)`) — by a **Dirichlet-form/capacity** trap. Every quantity is a min of the manifestly positive energy `∫|f′|²dν`; Cheeger `λ₁≥Φ²/4` supplies the gap from a **conductance floor**, never from a signed Hessian.

> **判词 (verdict).** The scheme is exact and it **does** dodge the indefinite-Hessian *operator*: no step evaluates `δ²K⪰δ` or signs `−D²+2adF_B`; the whole construction is positive energies and ratios of positive masses. But the floor uniform in (scale,`L`,cutoff,`k→0`) is **gap-strength**: on compact finite volume `Φ_*>0` is trivial (POTENTIAL-THEORY-CAPACITY §2), the content is uniformity, and at `k→0` the floor **is** H1 (`Var(a)≤1/λ`); the per-scale loss compounds and must be *summable* — the capacity face of BBS **marginal closure**. So: **operator bypassed; quantitative lock identical, key milder and better instrumented** (the free-field acid test is genuinely non-vacuous here).

## 1 · The object; why the Hessian never enters
Fix the coarse zero mode `a=⟨v,X⟩`, marginal `ν(da)=p(a)da`, `p∝e^{−V̄}`. Reversible Dirichlet form `E_ν(f)=∫|f′(a)|²dν`, generator `Lf=f″−V̄′f′`, spectral gap `λ₁(ν):=inf_{Var f>0}E_ν(f)/Var_ν(f)`. Cheeger + Poincaré:
`λ₁(ν) ≥ Φ_ν²/4`, `Φ_ν=inf_A Cap_ν(A,Aᶜ)/ν(A)`, `Cap_ν(A,Aᶜ)=inf{∫|f′|²dν : f≡1_A,0_{Aᶜ}}`;
`λ₁≥λ ⇒ Var_ν(a)≤E_ν(a)/λ=1/λ` (`E_ν(a)=1`) `⇒ G″(k→0)≳λ` (Ward `G″·Var(a)=1`, iff-Gaussian).
The energy `|f′|²⩾0`; the Savvidy depth `V̄″=−2g|F|` lives **in the measure**, never in the sign of the form — the quadratic form is `⩾0` for *any* `V̄`. Cheeger's RHS is (positive capacity)/(positive mass): **no indefinite operator is signed.**

## 2 · The trap
On the RG tower `j=0(IR)…J(UV)`, effective measures `ν_j` (block-spin marginals) and forms `K_j`. The old Hessian trap was `𝔗^H_j={‖R_j‖_X≤rg_j², δ²K_j⪰δ}` (FORGE5). The **C-BBS trap** is
`𝔗^C_j = { K_j : (i) ‖R_j‖_X≤rg_j² ; (ii) Dirichlet-energy density e_j(a):=|V̄_j′(a)|²≤M_j on the shell ; (iii) conductance floor Φ_{ν_j}≥Φ_*>0 }`.
This is a **ball in a positive norm `e_j=|force|²⩾0` plus a floor on a positive ratio** — no eigenvalue of an indefinite operator anywhere. What it controls: membership at every scale ⇒ the IR marginal `ν_0` inherits `Φ_0≥Φ_*` ⇒ `λ₁(ν_0)≥Φ_*²/4` ⇒ `Var(a)≤4/Φ_*²` ⇒ `G″≥δ`. The RG is not a sign map but the transport/propagation of a conductance floor down the tower.

## 3 · Four regions, what each controls
1. **Trap density (ii).** `e_j≤M_j` is the max-flow-dual certificate: it keeps `E(θ)=∫θ²/ν` finite so `Cap≥1/E(θ)` is constructive — the flow is expensive exactly where `ν` is thin. Sign-free accounting of the RG step.
2. **Small field** `|a|≲(N/12)^{−1/2}`. Prékopa (`Φ` convex, `Φ″⪰N/12`, PROOF-JACOBIAN-SAVES) ⇒ central marginal log-concave ⇒ core conductance `≥c(N/12)>0` (Cauchy–Schwarz/Brascamp–Lieb). The core cannot be the bottleneck.
3. **Large field, cone** `C_κ={2g|F|>μ²}`. **WSTC** (ABC-NEW-INEQ): `ν_j(C_κ)≤e^{−cβ(κ²−κ_typ²)N_□}`, uniformly in `L`. In flow language the cone is **measure-starved** and drops out of the relevant cut: the mass whose partition defines `Φ_ν` lives on the thermal shell, not the cone. *This is what neutralizes the large field* — an upper bound on mass, not a Hessian sign.
4. **Thermal shell** `|a|≤λ^{−1/2}`. The Savvidy dip is here, a saddle of `V̄` with measured `|V̄″|≈0.3`; barrier `H=|V̄″|·λ^{−1}≈0.3/λ=O(1)` at `λ=O(1)`. Eyring–Kramers ⇒ per-cut conductance `≥e^{−2H}` — a **finite** factor (≈0.55), no exponential suppression at the target scale. This is the only bottleneck, and it is `O(1)`.

## 4 · What remains to prove — exact residual
- **(R1) RG-invariance of the floor:** `T_j(𝔗^C_j)⊂𝔗^C_{j+1}` with `Φ_*` preserved. Coarse-graining can only *decrease* `λ₁`; the per-block loss is the Arrhenius factor `e^{−2H_j}`, so need `Σ_j H_j<∞` (`H_j=|V̄_j″|λ_j^{−1}`) — the **marginal log-resummation** of capacity loss, d=4 marginal ⇒ a *summability*, the capacity analogue of BBS marginal closure.
- **(R2) No thin neck:** certify the dip is a **bulk saddle** (`H` finite), not a narrow separating barrier with true conductance `→0`. Measured `0.3` supports; not proven.
- **(R3) Seed:** an a-priori `λ_0>0` at the UV (free-field / strong-coupling gap, unconditional: P3 `β<1/48`, SZZ) to start the bootstrap.
- **(R4) Uniformity:** `Φ_*` uniform in `L→∞` and in mode `→k=0`. This is where it becomes **H1**: `Var(a)≤4/Φ_*²` at `k→0`.

## 5 · Does it bypass the blocker? Honest
**Operator-level: yes.** No indefinite operator is signed; the FORGE5 obstruction (`remainder O(g²)` vs `cone depth O(μ²)`, no jet can flip the sign) vanishes because there is **no sign to flip** — the cone is killed by WSTC mass, the shell by a finite Arrhenius factor. **Quantitative-level: no.** (R1)–(R4) reproduce gap-strength content; (R4) at `k→0` is H1, and (R1) is the marginal problem. Gains are real but partial: (i) the **free-field acid test is non-vacuous** (flat zero mode ⇒ `Cap=0` ⇒ no bound), a genuine discriminator the Hessian route lacked; (ii) the input is an *asymptotic flow/mass* floor (strictly weaker than pointwise `G″≥λ`); (iii) the residual splits into four separately-attackable problems, the hardest (R1) being a **summability** rather than a **sign**.

**Sharpest line.** C-BBS trades "sign an indefinite Hessian on the cone" (impossible by scaling) for "lower-bound a convex Dirichlet energy" (signs nothing); the cone that sank FORGE5 is neutralized by its *own* measure-starvation, leaving a finite `O(1)` barrier on the thermal shell. But compounding that barrier per scale is the capacity face of marginal closure, and its `k→0` limit is H1 — **same lock, milder key, cleaner free-field handle.**

*Refs.* CONSTRUCTION-METHODS (c; blocker `G″_j=EH−VarG+W≥δ`); FORGE5-CONSTRUCTIVE (`δ²R_j=O(g²)`, cone `O(μ²)`); POTENTIAL-THEORY-CAPACITY (Wiener/equilibrium capacity, `λ₁≥Φ²/4`, finite-volume vacuity, Arrhenius `0.3/λ`); ABC-NEW-INEQ (WSTC cone starvation `e^{−cβΔκ²N_□}`); PROOF-JACOBIAN-SAVES (`Φ″⪰N/12`, Prékopa `W_geo≥0`, (◆)); G1-DRIVE / G2-DEPENDS-ON-G1 (L3 builds ν, H1 signs ν; `G″_a/a²≥Δ₀²`). External: Cheeger (1970); Jerrum–Sinclair, *SIAM J. Comput.* **18** (1989) 1149; Sinclair, *Combin. Probab. Comput.* **1** (1992) 351; Doyle–Snell, *Random Walks and Electric Networks*, MAA Carus **22** (1984); Eyring, *J. Chem. Phys.* **3** (1935) 107; Kramers, *Physica* **7** (1940) 284; Prékopa (1971); Bauerschmidt–Brydges–Slade, *LNM* 2242 (2019).
