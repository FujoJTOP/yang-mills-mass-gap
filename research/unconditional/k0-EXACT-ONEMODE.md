# k0-EXACT-ONEMODE — solve the constant mode exactly, then compare to 4D: the exact variance is real, the comparison inverts the sign

> **Verdict.** The constant mode **is** exactly solvable in tractable reductions — 2D YM (single-plaquette/holonomy, Migdal–Witten), the β=0 single link (Haar), the constant-field "cactus", the U(1) abelian-projected constant sector — each giving `Var(a₀)` as a closed finite-dim group integral, nonperturbative *within the model* and `L`-independent. But the **comparison does not close it: it fails on direction.** The only valid variance-transport is the total-variance/mixture identity `Var_full = E[Var(a|ξ)] + Var(E[a|ξ])`; the reduction freezes `ξ` and discards the nonnegative cross term. In the banked SCBL split this is exact: an exact one-mode model computes the **Gaussian–Schur curvature `S`**, so `Var_red = 1/S`, while the true `Var = 1/(S−D)`, `D≥0` the fiber Brascamp–Lieb defect. Hence `Var_red ≤ Var_true` and `gap_red ≥ gap_true`: the reduction yields an **upper** bound on the gap, exactly where a **lower** bound is required, and the deficit is `D = H1`. Not a closure.

## The exact objects (what genuinely closes)
- **2D YM, `T²`:** `ν(a₀)` = holonomy eigenvalue density, exact character sum; confining but **free-scaling** (`σ∝β`, no transmutation, no `a→0`).
- **β=0 single link (Haar):** `⟨tr U⟩=0`, `∫|tr U|²=1`, `⟨tr X²⟩` finite — `O(1)` variance *at the free point* (**spurious gap**).
- **Cactus:** `S_eff(a)=S_W(P(a))` at constant `a`; the *tree-level* curvature.
- **U(1) abelian projection:** exact constant-mode integral with monopole fugacity; Debye mass *if* monopoles condense.

Each is exact and nonperturbative-in-model; none is 4D SU(N) at β=2.5.

## The comparison — exactly why it inverts
Full marginal `ν(a)=∫ e^{−S_W(a,ξ)}J dξ`. Total-variance identity (exact):
```
Var_full(a₀) = E_ξ[Var(a₀|ξ)] + Var_ξ(E[a₀|ξ]).
```
The frozen reduction keeps the first piece at a point; the cross term `Var_ξ(E[a₀|ξ])` — the fiber-response of the constant mode — is **dropped** and is `≥0`. SCBL gives the same exactly: the one-mode model is the `D=0` (Gaussian) slice of `G″=S−D`; `Var_red=1/S`, `Var_true=1/(S−D)`. Since `D>0` in the Savvidy regime, `Var_true>Var_red`; the dropped term is the variance of `∂_ξ E[a₀|ξ]`. To reverse the inequality you must bound `D` above — **H1, verbatim**; a density-level comparison routes to a full-measure LSI = the gap (circular).

## Five-point verdict
1. **Nonperturbative:** ✓ *within* each reduction; ✗ as the target — all cases are 2D / β=0 / U(1) / cactus, none 4D SU(N) at 2.5.
2. **Sign-controlling:** ✗ — the reduction projects out the Savvidy-indefinite fiber that sets the 4D sign; model-sign does not control 4D curvature.
3. **Non-circular:** ✗ — valid transport needs `Var(E[a₀|ξ])≤·` = the force-variance/`D(v)` = H1.
4. **Non-SUSY + fails free field:** ✗ — Haar/2D/cactus give *finite* `Var` at the free point (spurious gap); they do not degenerate as the massless free field must (`Var(a_{k=0})=∞`, `S=0`).
5. **Uniform lower bound:** ✗ — by direction. Exact `Var_red` lower-bounds `Var_full`, so `gap_full ≤ gap_red=S` is an *upper* bound; a *lower* bound is required. The deficit `Var_true−Var_red` is exactly `D`.

## Exact failure
```
Var_red = 1/S ≤ 1/(S−D) = Var_true   (D = Var(f|0) − E⟨K^{−1}c,c⟩ ≥ 0),
⇒ gap_full ≤ gap_red = S.
```
The one-mode model pins `D=0`; the 4D gap is `S−D`; proving `S>D` is H1. Exact-one-mode + comparison therefore delivers an exact benchmark (`S`) and a **wrong-direction** bound — it cannot close the gap.

**Sharpest line.** The constant mode solves exactly, but every exact model *is* the Gaussian slice `D=0`; the only mix/marginalization inequalities run `Var_full ≥ Var_red`, so the exact variance lower-bounds the true one and the model gap upper-bounds the 4D gap — the sign the 5 tests forbid. The missing piece between the exact answer and the gap is precisely `D = Var(f|0)`, the same fiber Brascamp–Lieb defect that SCBL and PROOF-COARSE-POINCARÉ already isolate as H1. The direction of the only available comparison is the gap.

*Refs.* INVENT-VARIANCE-INEQ (`G″=S−D`, `D≥0`); PROOF-COARSE-POINCARÉ (`G″=E[H]−Var`, fiber indefinite); PROOF-JACOBIAN-SAVES (Prékopa, `W_geo`); PROOF-BORN-OPPENHEIMER (`V_eff″=G″`); PROOF-P3 (strong-coupling anchor, `β<1/48`); W7-Gpp0 (`G″≈1/Var`); W6-BE-at-2.5 (fiber negatives); Migdal, Sov. Phys. JETP **42** (1976) 413; Witten, CMP **141** (1991) 153 (2D YM exact solvability); Brascamp–Lieb 1976; Prékopa 1971; law of total variance; Strassen/convex order; Savvidy 1977.
