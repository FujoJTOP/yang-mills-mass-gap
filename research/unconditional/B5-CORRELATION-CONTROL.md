# B5-CORRELATION-CONTROL — SHAPE B: Elitzur decouples `a` from `∇S` only in the *pure-gauge* sector (tautological), and the Schmidt/Hessian block leaves the slow `k→0` block carrying both the gap and the Savvidy indefinite `−g|B|`; near-orthogonality is real but not a physical-sector theorem without H1

> **判词 (verdict).** B4/C4 fix the open-cell sign on the cross term `X=−2Cov(G,Φ′)` of `G″=Δ_dyn+W_geo+X`; the measured `corr(a,∇S)≈0` (with `corr(a,Φ′)≈+0.85`, so `Φ′≈c·a` by Stein) makes `X≈0`, `|X|≪W_geo`. Can this near-orthogonality be **proven structurally** — gauge symmetry (Elitzur) or a Schmidt/Hessian-block decomposition — to certify `|X|<W_geo` **without** a variance bound? **No: each structure is exact but silent exactly where the sign lives.** **(Elitzur.)** `S_W` is gauge-invariant ⇒ `∇S⊥` gauge orbits, and Elitzur kills gauge-variant expectations. Exact and nonperturbative — but only for the **pure-gauge component** of `a`, which is unphysical and must be projected out; the **transverse** sector (glueball gap *and* Savvidy instability) is untouched, and the measured `corr≈0` is largely this gauge-direction dilution. **(Schmidt.)** Splitting slow `k→0` ⊕ fast/stiff, the marginal curvature is the Schur complement `G″=A−B C^{-1}Bᵀ`; the stiff block is SZZ-gapped, but `A` — the slow-block Hessian — **is** `−D²+2adF`, whose `k→0` transverse eigenvalue is `−g|B|<0`; decoupling fast modes does **not** turn the slow block positive (that positivity is `W_geo`, Prékopa, already banked), and bounding the backreaction `BC^{-1}Bᵀ` (or inverting the indefinite `C`) is H1. Same fate as C5/C6, reached by another door. 5-pt: **1 PASS, 2 FAIL, 3 FAIL, 4 PASS, 5 FAIL.**

## §1 Elitzur route — decoupling in the wrong sector
`a=⟨v,X⟩`, `G=⟨v,∇S_W⟩`. Since `S_W` is a class function, `∇S_W` annihilates every gauge-orbit direction, and Elitzur forces `⟨F⟩=0` for gauge-variant `F`. Decompose `v=v_phys⊕v_gauge`: `G·v=G·v_phys` but `a=⟨v_phys,X⟩+⟨v_gauge,X⟩`. So the gauge part of the amplitude is *exactly* uncorrelated with the force, while the gap lives in `v_phys`. For the `k→0` soft mode the pure-gauge longitudinal component dominates `Var(a)`; the measured `corr(a,∇S)≈0` is **partly tautology**, and `corr(a,Φ′)≈0.85` confirms `Φ′` rides the gauge/large-`a` directions. Elitzur certifies orthogonality where nothing needs it and is silent on the deciding transverse sector.

## §2 Schmidt/Hessian-block route — the slow block is the wall
Integrate out fast modes: `G″=A−B C^{-1}Bᵀ` (Schur complement of `Hess V`). The fast/stiff `C` is gapped by block-positivity/SZZ (B3-ISOPERIMETRIC), so `C^{-1}` is controlled *there* — but `A` is the slow-block Hessian, and at `k→0` the transverse Savvidy eigenvalue is `−g|B|<0`. Momentum block-diagonality (`B=0`) is exact only for the free Gaussian field; the cubic `cos`-plaquette couples `k=0↔k≠0`, so `B≠0` nonperturbatively and `|BC^{-1}Bᵀ|<A` is the missing bound. This is verbatim B3-ISOPERIMETRIC's "`ρ_soft` *is* H1; cross-shell coupling unsigned."

## §3 Exact failing estimate
```
Elitzur:  G·v = G·v_phys ,  a = a_phys + a_gauge ,  Cov(a_gauge,G)=0 exactly
   ⇒ corr(a,∇S) small by dilution, NOT a physical-sector inequality; v_phys untouched.
Schmidt:  G″ = A − B C^{-1}Bᵀ ,  A@k→0 = −g|B| < 0 ;  C stiff = SZZ-gapped ✓
   ⇒ need |B C^{-1}Bᵀ| < A+W_geo ; B=0 Gaussian ⇒ free-field-exact (W_geo=0, no gap)
   ⇒ nonperturbative B bound, or C^{-1} on indefinite C  =  H1.
CS (open cell, β=2.5): |X| ≤ 2√(VarG·VarΦ′) ≈ 0.40 > W_geo ≈ 0.32 ⇒ G″ ≳ −0.08.
Closing |X|<W_geo ⟺ VarG<W_geo²/4VarΦ′ ≈ 0.8–1.3 (vs ≈2) = H1.   Ward: G″·Var(a)=1.
```

## §4 Five-point verdict (all five must hold)
1. **Nonperturbative — PASS (framework).** Elitzur + orbit-orthogonality are exact measure-level facts (no loops); the physical-sector control is the perturbative Schmidt decoupling, so this pass is largely nominal.
2. **Sign-controlling = LOWER bound — FAIL.** Neither route gives `G″≥λ` in the transverse sector: Elitzur is confined to `v_gauge`; `A` stays `−g|B|`.
3. **Non-circular — FAIL.** Bounding `BC^{-1}Bᵀ`, or inverting the indefinite `C`, is the Savvidy/H1 inverse; `VarG<1` ⟹ IR stiffness ⟹ H1 (Ward, zero slack).
4. **Non-SUSY + FAILS free field — PASS, but vacuous.** Momentum decoupling is **exact** for the free field, where `W_geo=0`, `corr=0`, `f₁²=0`; `|X|≤W_geo` returns `G″≥0` (no false positivity). But near-orthogonality is **free-field-saturated**, so it is not the YM sign-carrier — only `W_geo≥N/12` is YM-specific.
5. **Uniform `m_eff²(k→0)≥λ` — FAIL.** No uniform `λ`; the slow `k→0` Savvidy negativity is the target, undiminished by decoupling.

## Sharpest line
Gauge invariance and Schmidt decoupling both deliver **exact** orthogonality — in `v_gauge` (tautological, unphysical) and in the Gaussian momentum block (free-field-saturated). Neither touches the transverse `k→0` block, where `Hess V=−D²+2adF` still has `−g|B|<0` and `W_geo` (Prékopa) is the only positivity — already banked. The measured `corr(a,∇S)≈0` is real but **not a YM theorem**: a gauge-direction dilution plus a free-field Gaussian decoupling. Converting it into `|X|<W_geo` needs the slow-fast backreaction `BC^{-1}Bᵀ` bounded or the indefinite `C` inverted — verbatim H1. **Shape B closes the correlation door without moving the wall: the missing step is a nonperturbative, coupling-free bound on the transverse slow-mode backreaction — a species SU(N)+Wilson does not currently possess.**

*Refs.* B4-VARIANCE-MAG; C4-COV-LASTMILE (§5 tight `|X|`=H1); B3-TYPICAL-COERCIVITY; B3-ISOPERIMETRIC (stiff/soft split, `ρ_soft`=H1, cross-shell unsigned); C5-NEG-ASSOC; C6-DEPENDENCE-CLASS; PROOF-JACOBIAN-SAVES (`G″=Δ_dyn+W_geo+X`, `W_geo≥N/12`); PROOF-GAUSS-INVARIANT-PF; PROOF-SAVVIDY-EXACT-SIGN (`λ_min=−g|B|`); COMPUTE-GPP-SU3 (`corr(a,∇S)≈0`, `corr(a,Φ′)=0.85`, `E[Φ″]=0.458`, `W_geo=0.32`, `X=−0.04`); gauge-invariance/REPORT §3.2 (`∇S`⊥orbit, CCHS); glue-c-gauge-group §3.6 (Elitzur `F=χ(U_e)`→0). External: Elitzur, PRD **12** (1975) 3978; Osterwalder–Seiler, Ann. Phys. **110** (1978) 440; Shen–Zhu–Zhu, CMP **400** (2023); Brascamp–Lieb 1976; Prékopa 1971; Savvidy PLB **71** (1977) 133.
