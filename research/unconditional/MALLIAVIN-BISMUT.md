# MALLIAVIN-BISMUT — the stiffness `G″` as a Bismut–Elworthy–Li derivative of a conditional expectation: an exact identity whose Malliavin covariance *is* the defect `D`, entering with a minus sign — sign-defeating, reduces to H1

> **2026-11-23 · NEW SPECIES (Malliavin / Bismut–Elworthy–Li).** Forge: `G″(a)=d/da E[f|a]` is a *derivative of a conditional expectation*, so BE-L / the Malliavin integration-by-parts writes it as an expectation against a *conditional Malliavin weight* `W=−∂_a log q=f−m`, whose covariance `E[W²|a]=Var(f|a)` is the "stochastic-flow covariance". **Verdict.** The representation is *exact and nonperturbative* — and it is **sign-defeating in one line**: `G″(a)=E[∂_a f|a]+E[f·∂_a log q|a]=E[h_n|a]−E[W²|a]`. The flow covariance enters with a **minus**; the only positive Malliavin object is the *p-average* Fisher `∫G″p=J=Var(m)≥0`, which is H1. The closure it would need, `E[W²|0]≤E[h_n|0]−λ`, is the conditional Brascamp–Lieb of the fiber, valid iff `K=Π₀ Hess V Π₀⪰ρ>0` — exactly what Savvidy breaks at `k→0`. This is the *pathwise/Malliavin face* of D2 (Bakry–Émery, sharp via von Renesse–Sturm) and the *conditional* face of FORGE5-SL-MARTINGALE (QV = Malliavin covariance, residual `D`). **#1 PASS; #2,#3,#5 FAIL; #4 PASS (clean: free massless `G″=k²→0` exactly).**

## 1. The forge — conditional-Malliavin / Stein identity (exact, all β)
Chart `X=(a,y)`, `a=⟨v,X⟩`, fiber `y⊥v`; `μ∝e^{−V}`, `V=S_W+Φ`; `f=∂_a V`, `h_n=vᵀHess V v`; conditional density `q(a,y)=e^{−V(a,y)}/p(a)`; score `∂_a log q=∂_a log μ−∂_a log p=−f+m`, `m(a)=E[f|a]=G′(a)`. The conditional weight `W:=f−m=−∂_a log q` is `q`-centered. BE-L/Stein integration-by-parts gives
```
   G″(a) = ∂_a E[f|a] = E[∂_a f|a] + E[f·(∂_a log q)|a]
         = E[h_n|a] + E[f(m−f)|a] = E[h_n|a] − Var(f|a)
         = E[h_n|a] − E[W²|a].                                       (BE-L)
```
`E[W²|a]` is the **stochastic-flow / Malliavin covariance** of the force over the fiber. It is manifestly `≥0` and carries a **minus sign**. In pathwise form `Var(f|a)=∫E[vᵀA_t²v]dt` (`A_t`=fiber covariance process, `dA/dt`→Riccati) — literally the QV martingale covariance of FORGE5-SL-MARTINGALE. SCBL (`G″=S−D`) is the same identity after adding/subtracting the Schur term `⟨K⁻¹c,c⟩`: a *Gaussian–Schur* curvature `S` minus the Savvidy–BL defect `D`, `D>0 ⇔ K indefinite`.

## 2. Why it is sign-defeating, not sign-giving
BE-L outputs an **upper-capped** expression: `G″(a)=E[h_n|a]−E[W²|a] ≤ E[h_n|a]`. The flow covariance can only *lower* `G″`; it is the obstruction, not the tool. The one positivity the Malliavin calculus natively supplies is `E_p[G″]=∫(G′)²p = J = Var(m) ≥ 0` (Fisher, exact) — but this is the **p-average** of `G″`, not the point value `G″(0)`, and `J` is H1 in entropy clothing (E2). A "positively correlated flow" would need `E[W·(force-projection)]>0`; but `W` is the score, `E[W]=0`, and every pairing `E[f·W]` is a covariance of indefinite sign — only `E[W²]` is definite, and it appears subtracted.

## 3. Salvage — and its exact wall
The only way to turn (BE-L) into `G″≥λ` is a conditional Poincaré / Γ₂ on the fiber:
`E[W²|a]=Var(f|a) ≤ (1/ρ_fiber)E[|∇_y f|²|a]=(1/ρ_fiber)E[|H_{a⊥}|²|a]`, valid iff `K=Π₀ Hess V Π₀ ⪰ ρ_fiber>0`. Then `G″≥E[h_n|a]−(1/ρ_fiber)E[|H_{a⊥}|²|a]`. But at β=2.5 the fiber Hessian is **indefinite** (90–97% of Savvidy negatives in the fiber block, W6-BE), `ρ_fiber≤0`; the BL comparison *reverses* and injects `D≥0`. This is the "second-order Brascamp–Lieb transverse term" C4 §4 already found circular: the flow Jacobian `J_s` solves `dJ=−Hess V·J`, so its contraction rate is `inf spec(Hess V_eff)=G″` (von Renesse–Sturm, **sharp, zero slack**) — the Malliavin covariance *is* `Γ₂` *is* the Savvidy-indefinite Hessian.

## 4. Acid test — clean PASS
Free massless `k→0`: Gaussian `⇒ D=0` (Schur-tight, `Var(f|0)=0`), `S=1/Var(a)=k²→0`, so `G″=S−D=k²→0` **exactly** — BE-L manufactures *no* false gap (`f₁²=0`). Massive free: `G″=m²`. β=0 Haar: `Φ` convex `⇒ K⪰N/12 ⇒ D≤0`, `G″≥S≥N/12` — a theorem. Non-SUSY throughout. #4 PASS in letter *and* spirit (the free zero is genuine, not `Var`-tautological: `D=0` is the Gaussian statement, `S→0` the Schur statement).

## 5. Exact failing estimate
```
G″(0) = E[h_n|0] − E[W²|0]   (BE-L, EXACT; flow covariance subtracted)
  ⇒  G″(0) ≤ E[h_n|0]  (upper cap, wrong direction)
To close G″(0)≥λ need:   E[W²|0] = Var(f|0) ≤ E[h_n|0] − λ
  ⟺ conditional BL with  K = Π₀ Hess V Π₀ ⪰ ρ_fiber > 0  on the fiber
  ⟺ D(v) < S(v) − λ   [SCBL]   ⟺ Var(a_{k→0}) ≤ 1/λ   =  H1.
At β=2.5: K indefinite (90–97% fiber) ⇒ ρ_fiber≤0 ⇒ reversed BL ⇒ D>0; D/S≈2%, D<S measured, unprovable.
Only +Malliavin object: E_p[G″]=J=Var(m)≥0 (Fisher, p-average) → 0 free; = H1.
```

## 6. Five-point scorecard
1. **Nonperturbative — PASS.** Exact Stein/BE-L identity + pathwise QV, all β, no loop expansion.
2. **Sign-controlling `m_eff²(k→0)≥λ` — FAIL.** Flow covariance enters `−E[W²]≤0`; `E[h_n|0]` is Savvidy-indefinite; the rep gives `G″≤E[h_n]`, not `≥λ`.
3. **Non-circular — FAIL.** The only positive Malliavin quantity (Fisher average `J`, or nondegeneracy of `γ`) is H1 / gives nothing; the salvage needs `K⪰ρ` on the typical-action fiber = the Savvidy tachyon = H1.
4. **Non-SUSY + fails free field — PASS (cleanest).** Massless `G″=k²→0` exactly; massive `m²`; β=0 `≥N/12`.
5. **Uniform `≥λ` — FAIL.** `E[W²|0]≤E[h_n|0]−λ` is `D<S` is H1, nonuniform in `(L,a)`.

**Sharpest line.** Malliavin/BE-L *is* the natural language for `G″=d/da E[f|a]`, and it delivers the cleanest exact identity in the corpus — but the "stochastic-flow covariance" it adduces is precisely `Var(f|a)`, the defect `D`, and it arrives with a **minus**: BE-L caps `G″` from above, it cannot floor it. The species is sign-defeating by construction, and its salvage (conditional Brascamp–Lieb on `K=Π₀ Hess V Π₀`) is the Savvidy-indefinite fiber = H1. It therefore joins D2 (`λ=G″`, von Renesse–Sturm sharp) and FORGE5-SL-MARTINGALE (QV residual `=D`) as the third face of one object: **positivity needs `K⪰ρ>0` at typical action and high momentum, and that inequality is the gap.** Not 5/5.

*Refs.* Bismut, *Large Deviations and the Malliavin Calculus* (1984); Elworthy–Li, *Diff. Geom. Appl.* **4** (1994) 13 (BE-L); Malliavin (1978); Nualart, *The Malliavin Calculus* (2006); Fournier (2002); Bakry–Émery 1985; von Renesse–Sturm, CPAM **58** (2005) 923. Internal: FORGE5-SL-MARTINGALE (QV `Var(a)=∫E[vᵀA_t²v]`, residual `D`); D2-COUPLING-DOEBLIN (BE, `λ=G″` sharp); INVENT-VARIANCE-INEQ (SCBL `G″=S−D`); ABC-NEW-INEQ (`G″=E[H]−Var`, Ward `E[af]=1`); OTHER-H1-POINCARE (`Hess V_eff=G″=E[PᵀHP]−Var`, two-route exhaustion); PROOF-STOCHASTIC-LOCALIZATION (`A_t=(H+tI)⁻¹`); PROOF-JACOBIAN-SAVES; W6-BE-at-2.5 (fiber indefinite, 90–97%); C4-COV-LASTMILE (§4 second-order BL circular); CONTINUUM-VARIANCE-DEFECT (`D/S` IR denominator = H1).
