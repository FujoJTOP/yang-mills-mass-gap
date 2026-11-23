# D-BBS-MEASURE-TRAP — replace "convex at large field" by "measure-starved at large field": the two-region trap is flow-invariant and acid-clean, but the thermal shell keeps an O(1) Arrhenius barrier that only becomes H1 at k→0

> 2026-11-23 · constructive. The BBS large-field leg works for φ⁴ because the large field is *more* convex (`+12λφ²` makes the LF Hessian positive). For gauge the large field is *less* convex: the Savvidy Hessian `−D²+2 ad F_B` is indefinite on the cone. **New mechanism:** do not ask the large field to be convex — ask it to be *absent*. Replace BBS's LF Hessian positivity by the proven **Ward–Stein transport-cone (WSTC) measure starvation** `ν(C_κ)≤e^{−cβ(κ²−κ_typ²)N_□}`. Build the two-region trap (small field log-concave ∪ large field starved), test flow-invariance (the BBS requirement), and quantify the uncovered transition shell.
>
> **判词 (verdict).** The two-region trap is **real and flow-invariant** — and in a *stronger* sense than BBS: each region-leg is re-proven at every scale (Prékopa's `Φ″⪰N/12` is scale-free; WSTC is per-scale and *strengthens* into the UV), so `T_j(𝔗_j)⊂𝔗_{j+1}` holds by re-establishment, not contraction. Neither leg signs an indefinite operator — the core signs a convex form, the cone is killed by *mass*. This genuinely ports BBS's large-field leg to gauge. **But trap-invariance is not gap-invariance:** the legs cover core and cone and leave the **thermal shell**. At one scale that shell is a **finite Arrhenius `O(1)` barrier** (`H≈0.3/λ`, per-cut conductance `≈0.55`), *not* H1. H1 survives only in the **k→0 uniformization** of that barrier (no-thin-neck `R2` + `Var(a)≤1/λ` `R4`). The cone is neutralized; the gap is one uniformization short.

## 1 · The two-region trap (a)
Partition field space at scale `j` into `Ω = Core ⊔ Shell ⊔ Cone` and define
```
𝔗^MS_j = { ν_j :  (C-core)  Hess(S_W+Φ)_vv ⪰ c(μ_j²+N/12) > 0 on {‖F‖≤a_c},
                    (C-starve) ν_j(C_κ) ≤ e^{−cβ_j(κ²−κ_typ²)N_□},  C_κ={‖F‖²≥κ²}>Cone,
                    (C-balance) ν_j(Shell) ≤ η_j }.
```
- **Core** `{‖F‖≤a_c}`: `−D²⪰μ_j²`, Savvidy term small, and `Φ″⪰N/12` (SU2 `1/6`, SU3 `1/4`; PROOF-JACOBIAN-SAVES). So `Hess⪰c(μ_j²+N/12)>0` ⇒ marginal log-concave ⇒ Brascamp–Lieb floor `Var(f)≤E|∇f|²/curv`. **真定理**, scale-free via Prékopa.
- **Cone** `{‖F‖²≥κ²}`: **WSTC** ⇒ mass `≤e^{−cβ_j(κ²−κ_typ²)N_□}`, uniformly in `L`. **真定理** (ABC-NEW-INEQ §2).
- Neither leg evaluates `δ²K_j⪰δ` or signs `−D²+2ad F_B`: one is a convex curvature, the other an *upper* bound on measure.

## 2 · Is the trap flow-invariant? (b)
BBS requires `T_j(𝔗_j)⊂𝔗_{j+1}`. Here both legs are *re-established at every scale*, which is stronger than contraction:
- `(C-core)_{j+1}`: as `j↑` (UV), `g_j²↓` and `μ_j²↑`, so the foreground floor **rises**; `Φ″⪰N/12` is pointwise and scale-independent. Holds irrespective of `T_j`.
- `(C-starve)_{j+1}`: WSTC is proved from `E_μ[fa]=1` + large deviations, no gap input; `β_j=2N/g_j²↑` into the UV ⇒ the exponent **strengthens**. Holds irrespective of `T_j`.
So `𝔗^MS` is invariant, indeed monotone-strengthening. This is *cleaner* than BBS, which needs LF convexity to contract **both** halves; gauge needs it on neither — the cone half is removed by mass, not curvature.
**Honest catch:** trap-invariance `≠` gap-invariance. The gap is the *global* conductance floor `Φ_ν=inf_A Cap(A,A^c)/ν(A)` (`λ₁≥Φ²/4`). The legs fix cuts inside the core and into the starved cone, but the **bottleneck cut is a balanced one in the uncovered Shell**. Invariance of the pair is necessary, not sufficient.

## 3 · The residual: the thermal shell (c)
The Shell is the transition `|a|~λ^{−1/2}` where **neither** estimate applies: log-concavity fails (Savvidy dip `V̄″=−2ad F_B<0`, `|V̄″|≈0.3`), and starvation fails (the Shell carries the `O(1)` bulk mass, not `e^{−c…}`). It is the only uncovered region. Because the dip is a *saddle* of `V̄`, crossing it costs the Eyring–Kramers barrier
```
H = |V̄″|·λ^{−1} ≈ 0.3/λ = O(1)   at   λ=O(1),   per-cut conductance ≥ e^{−2H} ≈ 0.55.
```
**Verdict on (c): at one scale the Shell residual is the Arrhenius `O(1)` barrier, NOT the full H1.** That is the exact gain — measure-starving the cone leaves only a *finite* thermal barrier to cross, replacing FORGE5's impossible "sign an indefinite Hessian on the cone." H1 lives in three caveats:
- **(R1) Compounding.** The per-block loss `e^{−2H_j}` multiplies over `J~log(1/a)` scales. With `H_j~g_j²` this *is* the marginal log; by C-BBS-MARGINALLOG (A) it resums into the running coupling given L2, and the deviation `Σg_j⁴<∞`. **Ports — not a new wall.**
- **(R2) No thin neck.** `0.3` is measured; the dip must be a *bulk* saddle (`H` finite), not a narrow separating barrier with true conductance `→0`. Unproven.
- **(R4) k→0 uniformity.** At the softest mode the Shell and the whole volume coincide; the Arrhenius estimate presupposes a width/distribution it cannot supply, and the floor collapses to `Var(a)≤1/λ` — **exactly H1**.

## 4 · Ledger and sharpest line
| leg | status |
|---|---|
| (C-core) log-concave floor `N/12`, scale-free | **真定理** (Prékopa) |
| (C-starve) cone `e^{−cβΔκ²N_□}`, per-scale, `L`-uniform | **真定理** (WSTC) |
| `T_j`-invariance of `𝔗^MS` (re-establishment) | **真定理** |
| thermal-shell Arrhenius `H≈0.3/λ=O(1)` | **条件定理** (needs R2) |
| `ΣH_j` compounding | **ports given L2** (marginal log) |
| k→0 uniformization of the barrier | **OPEN = H1** |

**Acid test.** Free massless field: `F=0` ⇒ cone empty (WSTC→0), no Haar ⇒ `Φ″=0` ⇒ no core floor, `Var(a)=∞` ⇒ no gap. The trap is **non-vacuous**: it bounds YM, not the free field.

**Sharpest line.** Replace "convex at large field" with "starved at large field" and the Savvidy cone stops being the obstacle — a **true theorem kills it by mass, not curvature** — leaving only a finite `O(1)` thermal barrier. The trap is invariant; H1 was relocated, not removed, to the `k→0` uniformization of that barrier. The cone is the wrong enemy; the thermal shell's dilute limit is the right one.

*Refs.* ABC-NEW-INEQ (WSTC cone starvation); C-BBS-FORMULATE (capacity trap, Arrhenius `0.3/λ`, `E(θ)=∫θ²/ν`); C-BBS-MARGINALLOG §1–4 (marginal-log ports, LF convexity breaks); C-BBS-UVTRAP (per-shell floor rises into UV); FORGE5-CONSTRUCTIVE (`δ²R_j=O(g²)` vs cone `O(μ²)`); PROOF-JACOBIAN-SAVES (`Φ″⪰N/12`, Prékopa `W_geo≥0`); CONSTRUCTION-METHODS (c; blocker `G″_j=EH−VarG+W≥δ`); PROOF-TRAPPING-REGION §1. External: Bauerschmidt–Brydges–Slade, *LNM* 2242 (2019); Eyring, *J. Chem. Phys.* **3** (1935) 107; Kramers, *Physica* **7** (1940) 284; Prékopa (1971); Brascamp–Lieb (1976); Savvidy, *Phys. Lett. B* **71** (1977) 133.
