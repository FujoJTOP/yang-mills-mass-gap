# H2-GLOBALIZATION — block-local ⇒ full-fiber Poincaré/LSI, without a contraction

**Verdict.** Globalization to the full fiber at β=2.5 **does not close**, and the four proposed routes die on one and the same estimate — the **inter-block (renormalized) coupling** — *not* the intra-block Savvidy UV concavity. A full-fiber Poincaré/LSI **does** hold at every finite volume, but its constant is `e^{2βN_□(V)}C_G` — **not uniform in V**. So H2 in the *full-fiber* form CONDITIONAL-LATTICE-GAP Step 2 uses is **open**; only the block-local H2 (PROVE-H2, `Λ_UV=e^{−2βN_p(B)}/C_G`) is proved.

## Why the local gap does not tensorize
The fiber ξ=(ξ₁,…,ξ_{n_B}) carries `μ(·|U′)∝∏_block Haar·e^{−βH_fib}`, and H_fib couples neighbouring blocks through the plaquettes **straddling their shared faces**. Were the fiber the product `⊗μ(ξ_i|U′)`, Poincaré would **tensorize** with no loss:
`Var_μ(F) ≤ (1/Λ_UV) Σ_i E_μ|∇_{⊥,i}F|²`.
The blocks are not independent: the inter-block oscillation is `2β·(#shared plaquettes)=O(β)≫1`. Tensorization needs zero coupling; Holley–Stroock comparison needs *small* coupling; Zegarlinski localization needs a *contraction*. None holds at β=2.5.

## (A) Trivial HS globalization — rigorous, not V-uniform
Against `⊗Haar` over blocks, `osc_ξ(βH_fib) ≤ 2βN_□(fiber)`, `N_□(fiber)=Σ_blocks N_p(B)~V`. Holley–Stroock:
`Ent_μ(f²) ≤ e^{2βN_□(fiber)}C_G·E|∇_⊥f|²`, and the same for Var.
This is the only *rigorous* full-fiber inequality available here; its constant is `e^{cV}`-bad. **FAILS uniformity.**

## (B) Zegarlinski / Nachtergaele localization — uniform iff the contraction/LGP holds
Let `a_{ij}=osc_{ξ_j}(log π^{cond}_i)`, `C(β)=sup_i Σ_{j~i} a_{ij}`. The one-site-to-global (Zegarlinski) and martingale (Lu–Yau; Nachtergaele's concentric-block) theorems both give a volume-uniform gap
`σ_glob ≤ σ_UV/(1−C(β))` (Zegarlinski), or `Δ≥γ` when the per-scale conditional gap `γ_j≥γ` uniformly (Nachtergaele),
**only if** `C(β)<1` / `γ_j≥γ` uniformly in scale and boundary condition. At β=2.5 the raw log-oscillation is `2β=5` per shared plaquette, so `a_{ij}~2β·(#shared)` and `C(β)≫1` (a single coupling already gives `C>1`) ⇒ `1−C<0`, the Neumann series diverges and the bound grows with the torus cycle length. Equivalently, PROVE-H2 supplies Nachtergaele's `γ_j` **only at scale 0**; for growing blocks `γ_j` is the boundary-uniform conditional gap = LGP = the k→0 transverse mass (IDEA-finite-size), which VAR-FINAL leaves uncontrolled (`G″(0)≈−0.165`). **FAILS at the block coupling.**

## The four methods are one wall
- **Martingale / Lu–Yau and (3) conditional telescoping.** Peeling blocks in any order, `μ(ξ_k|ξ_{<k})` is a *mixture* of block-conditionals over the integrated-out future; mixtures destroy the Poincaré constant (two well-separated bumps: Var∼L², gradient budget O(1)). The renormalized marginal is uncontrolled — block coupling.
- **(2) Multiscale Dobrushin.** One b=2a step leaves `β_eff≈β≫1`; reaching `β_eff<1` is the deep-IR RG flow (IN-9). Block coupling.
- **(4) RP + additivity.** RP (OS'78) is a time-reflection sign; plaquettes straddle block faces, so the fiber generator is not a sum of commuting block generators and gaps do not add (P4). Block coupling.

## The failing estimate, pinned
**It is the block coupling, not Savvidy.** Savvidy concavity is *intra-block*, already absorbed by each block's bounded-oscillation HS gap; an HS-type localization is blind to Hessian sign. Globalization needs the *inter-block influence sum* `C(β)<1` (equivalently per-scale `γ_j≥γ`), which β=2.5 violates. Savvidy obstructs only the *alternative* global-convexity/BE route (W6) — a different road.

## Consequence — and why it is IR, not technical
CONDITIONAL Step 2 bounds the **whole-fiber** `Var_{μ(·|U′)}` by `(1/Λ_UV)E|∇_⊥|²`; block-local H2 gives this only under tensorization, false here, so Step 2 requires exactly this open globalization. A contraction-free globalization would grant volume-uniform gaps to *all* bounded-oscillation cyclic Gibbs systems, contradicting the strong-coupling-or-structure necessity (P4): the gap between block-local and full-fiber **is** the IR/mass-gap content. Volume-uniform: **NO**.

*Refs.* Holley–Stroock 1987; Zegarlinski CMP 1989; Guionnet–Zegarlinski; Lu–Yau 1993; Nachtergaele CMP 175 (1996); Nachtergaele–Sims MPRF 12 (2006); OS'78; W6-BE; W4; VAR-FINAL; IDEA-finite-size; PROVE-H2; P4.
