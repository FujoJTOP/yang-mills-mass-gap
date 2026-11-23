# K0-OVERLAP-AUDIT — is the k=0 coarse mode in the Savvidy cone, or is the −0.3 a different (vector) mode a 0⁺⁺ projection decouples?

> 2026-11-23 · PRECISE OVERLAP AUDIT of the NUMERICS-AUDIT tension: conditioned `G″=+0.4` vs bare `PᵀHP` min `−0.3`. **判词 (verdict).** The "coarse mode" the campaign actually measures is a **single transverse plane-wave gluon amplitude** `a=⟨v,log U⟩` (one axis × one polarisation × one colour; `k0_continuum_run.py` `v=plane_mode(L,ax,pl,c,n)`) — gauge-*dependent*, colour-*adjoint*, a `1⁻` vector. **It is not the gauge-invariant 0⁺⁺ scalar.** So (a) it sits in the **same** transverse-gluon sector as the Savvidy negatives — overlap is *not* quantum-number-suppressed; (b) the literal 0⁺⁺ scalar is genuinely orthogonal to the negative vector modes, but (c) *no campaign number measures it*: the `+0.3` and the `−0.3` are the **conditioned vs bare versions of one gluon**, not a 0⁺⁺-vs-vector pair. The negative is **not** "removed by conditioning": conditioning *adds* `W_geo`, and the sign flip is an **uncertified value cancellation** (LB `−0.10`). The 0⁺⁺ decoupling is real but does **not** discharge H1, which is measure-level (adjoint) and *includes* the cone.

## (0) Three objects — the numbers do not measure their labels

| label | object | sector / QN | gauge | what it is |
|---|---|---|---|---|
| `+0.307`/`+0.192` (cond. `G″`) | `a=⟨v,log U⟩`, `v`=plane wave `(ax,pol,c,n)` | single transverse gluon, adjoint, `1⁻` | **dependent** | curvature of `−log P(a)`, `L=8`, `k=π/4` |
| `−0.3` (`PᵀHP` min) | coarse-restricted **bare** Wilson Hessian | transverse-gluon, low-momentum block | covariant (Hessian of `S_W`) | bare `E[H|coarse]` min, β=2.5 |
| 0⁺⁺ "glueball" | plaquette field `Σ ½Re tr U_p` (`corr_su2.py`) | color-**singlet** scalar `0⁺⁺`, bilinear | **invariant** | 2-pt correlator → `ξ` ("estimate, not proof") |

The first two are the **same operator sector**; the third is the physical 0⁺⁺ and is a *different* object measured only as a correlation length, never as a curvature sign.

## (a) Overlap ε — the coarse mode is a gluon, and it overlaps the cone
Savvidy negatives are transverse 1-forms `δA_⊥^a` (adjoint, vector), `n=0` Landau level, `s=−` magnetic, `k_∥=0` — the **same** adjoint/vector sector as the measured transverse plane wave of `log U`, so energy-orthogonality is no defence. Two overlap facts, both true:
- **Analytic (pure constant-`B`): `ε=O(1)`.** The band `λ_{0,−}=k_∥²−g|B|` bottoms at `k_∥=0`; a `k_⊥=0` plane wave has large weight on the width-`ℓ=1/√(gB)` Landau orbital (D-BBS §c).
- **Empirical (β=2.5): small for the population.** The most-negative *fine*-Hessian eigenvector is `90–97%` in the UV/fluctuation block, only `3–10%` in the coarse collective block (W6-BE §2). Yet a **softened coarse remnant survives blocking**: `PᵀHP` min `≈−0.3` (PROOF-COARSE-POINCARE §a).

So ε is neither zero nor `O(1)` for the bulk — `O(0.03–0.1)` in population, plus a **residual coarse negative of depth `−0.3`**. The k=0 coarse gluon **is in the cone**, lightly.

## (b) Quantum numbers — the 0⁺⁺ *is* orthogonal, but it is not what was measured
`⟨0⁺⁺|single transverse gluon⟩=0` by colour (singlet vs adjoint), `C`-parity (`+` vs `1⁻`), and composite-vs-elementary — genuine decoupling. **But the campaign's `+0.3` is a gluon, not the 0⁺⁺.** The clean "scalar positive / vector negative" split is **not what the data compare** — both numbers are vector.

## (c) Which mode owns `−0.3`, and why "project onto 0⁺⁺" does not rescue H1
The `−0.3` is the **bare coarse-restricted transverse-gluon** Hessian, same sector as the `+0.3` it competes with. Apples-to-apples at β=2.5: bare fine `min=−2.8`, coarse `−0.3`; exact conditioned `G″=+0.51` (PROOF-JACOBIAN-SAVES §3). The `−0.3→+0.3` flip is **not** a projection (the mode never changes sector); it is `W_geo=E[Φ″]−Var(Φ′)` beating `|E[H]|+Var`. That term is proved only `≥0` (Prékopa), value-free; keeping the certified piece gives LB `−0.10` — **0 not excluded**.

The tempting escape — "−0.3 is a vector mode, so A(k=0) on the 0⁺⁺ scalar is positive" — is **true as an operator statement and useless as a proof**:
1. The campaign's "`G″(k→0)>0`" evidence is the **gluon** (`+0.3`), same sector as the cone — there is no 0⁺⁺ curvature measurement in COMPUTE-GPP/NUMERICS-AUDIT.
2. H1 (PROOF-IR-H1) is a **measure-level coarse Poincaré on `∇_∥`** — an LSI of the orbit-space measure controlling *all* gauge-invariant observables including 0⁺⁺. An LSI cannot be restricted to the 0⁺⁺ subspace: the negative transverse directions degrade the measure's log-concavity, and the inequality needs the whole function space.
3. "A(k=0) on the 0⁺⁺ positive" is the **glueball mass from a gauge-invariant 2-pt correlator** (`corr_su2.py`): lattice `~1.7 GeV`, well-known **evidence**, but a *different* computation (decay rate, not curvature sign) that never feeds the SLC/H1 route.

## Sharpest line
The measured coarse mode is a **transverse gluon**, not the 0⁺⁺ scalar; it overlaps the cone at `3–10%` of the negative population with a **residual bare `−0.3`** — so it is **in the cone**, not decoupled. The physical 0⁺⁺ **is** orthogonal to the negative vectors, but that concerns a two-point operator **no campaign run measured**, and it does **not** remove the measure-level obstruction: conditioning does not *project out* the negative, it *amortises* it with an uncertified Haar term — a cancellation, not a decoupling (LB `−0.10`). Honest: the 0⁺⁺ scalar is **decoupled**; the coarse mode carrying `G″`/H1 is **not**, and **H1 is the gap**.

*Refs.* `k0_continuum_run.py` (`v=plane_mode`, `a=⟨v,log U⟩`); `corr_su2.py` (gauge-invariant plaquette 0⁺⁺-like `ξ`); W6-BE-at-2.5 (§2 negatives `90–97%` UV/fiber); D-BBS-LOWRANK-KATO (§c `k_∥=0` soft, `ε=O(1)`); PROOF-COARSE-POINCARE (§a `PᵀHP` min `−0.3`); PROOF-IR-H1 (H1 on `∇_∥`, rename); PROOF-JACOBIAN-SAVES (§3 three objects, `W_geo≥0`, LB); PROOF-SAVVIDY-EXACT-SIGN (§3: negatives are *not* "ASD" — SD/ASD both BPS/PSD, instability needs non-extremality, max at `Q=0`); NUMERICS-AUDIT; COMPUTE-GPP-K0-CONTINUUM (LB `−0.098`).
