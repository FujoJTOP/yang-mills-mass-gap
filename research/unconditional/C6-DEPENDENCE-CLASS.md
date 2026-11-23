# C6-DEPENDENCE-CLASS — SHAPE C round 5: the dependence-inequality landscape is a strict dichotomy; every right-sign class collapses to Cauchy–Schwarz, every survivor is sign-blind or wrong-signed, and none passes the 5

> **判词 (verdict).** We need either a **negative** cross-sign `Cov(G,Φ′|0)≤0` (⇒ `X=−2Cov≥0` ⇒ `G″≥W_geo+Δ_dyn≥N/12`, no variance) **or** a variance collapse `VarG≤cW_geo`, for SU(N) Haar×Wilson on the open cell (β=2.5, typical field, high-momentum). The landscape splits cleanly, and the split *is* the verdict. **(Bin 1 — right sign, unsatisfiable.)** NA, NR/NU/NOD, strong Rayleigh, log-submodular, DPP/Pfaffian, Chebyshev/TP all *aim* at `Cov≤0`, but each presupposes an order, a real-stable polynomial, a DPP kernel, disjoint monotone coordinates, or a total order on the `2(N²−1)V`-dim fiber that compact SU(N) on a delocalized slice lacks; the only one mature enough to degenerate, Chebyshev, degenerates *exactly* to Cauchy–Schwarz (`|X|≤2√(VarG·VarΦ′)≈0.40 > W_geo≈0.32`). **(Bin 2 — magnitude, circular or wrong-signed.)** FKG/Griffiths is the **opposite** sign (`Cov≥0⇒X≤0`, anti-helpful) and also needs an absent order; Brascamp–Lieb needs `Hess V⪰0` — the Savvidy wall itself (`λ_min=−g|B|<0`) — and returns `VarG≤cW_geo`=H1; Simon–Lieb is a sign-blind variance bound, =H1 on the coarse slice. **Only Prékopa (`W_geo≥0`) holds on the cell — and that is the banked `N/12` floor, not the cross-sign `X`.** Nothing passes #2/#3/#5.

## Structural locus of all dependence
In the flat chart the Haar base is a **product over links**: `dU_ℓ=J(X_ℓ)dX_ℓ`, so `Φ=−log J=Σ_ℓ φ(X_ℓ)` is a single-link sum and the β=0 measure is `⊗_ℓ e^{−φ(X_ℓ)}dX_ℓ` (independent links, `Var_Φ(a)=O(1)`). Hence **all inter-link dependence — the only source of `Cov(G,Φ′|0)≠0` beyond local same-link terms — lives in the Wilson factor `e^{−S_W}`**: a class function of `tr U_□`, flat on gauge orbits, with no order-compatible (super/sub)modularity on SU(N) (compact ⇒ torsion ⇒ no lattice-ordered group; PROOF-COMPARISON-FKG §1). Every association class is thus a question about `e^{−S_W}`.

## Per-candidate table
| Class | Holds on cell? | Sign it gives | Target | Load-bearing absent hypothesis |
|---|---|---|---|---|
| Negative association (NA) | **No** | ≤0 (right) | `Cov≤0` | no lattice order; `G,Φ′` gauge-covariant, delocalized, overlapping — not disjoint increasing coordinates; `a=0` breaks product |
| Neg. dependence NR/NU/NOD | **No** | ≤0 (right) | `Cov≤0` | needs orthants/order; forces are smooth+delocalized, not indicator events |
| Strong Rayleigh | **No** | ≤0 (right) | `Cov≤0` | no real-stable generating polynomial over continuous group configs |
| Log-submodularity | **No** | ≤0 (right) | `Cov≤0` | no lattice; `e^{−S_W}` flat on orbits ⇒ not sub/super-modular |
| DPP / Pfaffian / PPP | **No** | ≤0 (DPP) but physics **≥0** | `Cov≤0` | no kernel rep; Savvidy fiber is *attractive* (`D>0`), the reverse of repulsion |
| Brascamp–Lieb (2nd-order) | **No** at open cell | magnitude only | `VarG≤cW_geo` | needs `Hess V⪰0` = Savvidy wall (indefinite); else optimal `B=Hess⁻¹`=H1 |
| FKG / Griffiths | **No** | **≥0 ⇒ X≤0 (WRONG)** | anti-target | no order **and** wrong sign: lowers `G″` |
| Simon–Lieb | vacuous/circular | none (sign-blind) | `VarG≤cW_geo` | variance magnitude on coarse slice = susceptibility = H1 |
| Chebyshev / total-positivity | **No** | ≤0 (right) | `Cov≤0`→CS | no total order / common scalar sufficient statistic; multivar rearrangement=Fréchet–Hoeffding=CS |

## Five-point scorecard (covers every class)
1. **Nonperturbative — PASS** (measure-level; no loops/`1/N`/`g`). 4. **Non-SUSY + fails free field — PASS** (`F≡0⇒G=Φ′=0⇒Cov=0`, no spurious positivity).
2. **Sign-controlling lower bound — FAIL** (Bin 1 reduces to CS, `G″≥−0.08<0`; Bin 2 sign-blind/wrong-signed). 3. **Non-circular — FAIL** (variance output = `Var(a)≤1/λ`=H1; BL's optimal `B`=indefinite inverse; `Cov≤0` on-cell needs the copula=measure). 5. **Uniform `m_eff²(k→0)≥λ` — FAIL** (no `λ`; Bin 1 no sign theorem, Bin 2's constant is H1).

## Exact failing estimate
```
Need  X ≥ λ−W_geo ≈ −0.22  (⟺ Cov(G,Φ′|0) ≤ 0.11)   or   VarG ≤ W_geo²/(4VarΦ′) ≈ 0.8–1.3.
Bin 1 (NA/Chebyshev/log-submod./DPP): aims Cov≤0 — collapses to
   Fréchet–Hoeffding  |X| ≤ 2√(VarG·VarΦ′) ≈ 0.40 > W_geo ≈ 0.32  ⇒ G″ ≥ −0.08 < 0.
Bin 2 (BL/Simon–Lieb): VarG≤c  *is*  Var(a)≤1/λ = H1; BL needs Hess V⪰0 (false: λ_min=−g|B|).
FKG/Griffiths: Cov≥0 ⇒ X≤0 ⇒ G″ ≤ W_geo+Δ_dyn  — WRONG direction.
Per-mode:  X∈[−0.09,+0.16], mean≈0  ⇒ no uniform negative covariance exists to certify.
```

**Sharpest line.** The landscape splits, and the split is the answer. Every class that would deliver `Cov≤0` needs a structure SU(N) lacks — an order (torsion forbids a lattice-ordered group), a real-stable polynomial, a DPP kernel, or a total order on a `2(N²−1)V`-dim fiber — and the one mature enough to degenerate, Chebyshev/NQD, degenerates *exactly* to Cauchy–Schwarz, the 0.40-vs-0.32 shortfall already banked. Every class that survives without an order is either wrong-signed (FKG/Griffiths push `Cov` up) or sign-blind and equal to H1 (BL/Simon–Lieb variance = the susceptibility ceiling). The physical anti-localization (`G∝|F|` at O(π) links; `Φ′` surge at the `F≈0` log-cut) is a *mechanism*, not an association theorem — and per-mode data (`X` sign-varying) show it does not uniformize. **No class closes the last mile; the required object is a deterministic, coupling-free negative-correlation comparison for `e^{−S_W}` — a species that does not yet exist.**

*Refs.* PROOF-COMPARISON-FKG (§1 no order, §2 wrong sign); C4-COV-LASTMILE (§2–§5); B3-TYPICAL-COERCIVITY; B3-FUNCTIONAL-BE; ABC2-CONSTRUCTIVE §4; INVENT-VARIANCE-INEQ (BL needs `K⪰0`, `D>0` attractive); PROOF-SAVVIDY-WGEO-DOMINATION; PROOF-JACOBIAN-SAVES; OTHER-SAVVIDY-SIGN §2 (product Haar); COMPUTE-GPP-SU3. External: Harris 1960; Griffiths 1967; Fortuin–Kasteleyn–Ginibre 1971; Pemantle, J. Math. Phys. **41** (2000) 1371; Borcea–Brändén–Liggett, JAMS **22** (2009) 521; Lyons, Publ. Math. IHÉS **98** (2003) 167; Brascamp–Lieb 1976; Prékopa 1971; Fréchet 1951/Hoeffding 1940/Rüschendorf 1983; Karlin, *Total Positivity* (1968); Simon–Lieb (conditional-variance decomposition; cf. Efron–Stein 1981); Savvidy, PLB **71** (1977) 133.
