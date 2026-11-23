# RCD-SYNTHETIC-RICCI — the gauge measure space `(A/G, g, μ)` as an `RCD(K,∞)` space: synthetic Ricci is **not weaker** than pointwise `Hess⪰K` (on flat `A` they are Bakry–Émery-equivalent), the one positive term (O'Neill `|F|²`) averages to the gluon condensate, and the only positive Poincaré is the compactness artifact that dies as `L→∞`

> **2026-11-23 · NEW SPECIES — synthetic Ricci curvature / Lott–Villani–Sturm `RCD(K,∞)`.** 判词 (verdict). View `(A/G, g, μ)` (`A=su(N)^V` flat, `μ=Haar·e^{−S_W}`, orbital quotient, coarse amplitude `a=⟨v,X⟩`) as a metric-measure space and seek `RCD(K,∞)` with `K>0`: it yields `Var(f_t)≤e^{−2Kt}Var(f)` ⇒ Poincaré `λ₁≥K` (Savaré 2014) ⇒ `Var(a)≤(1/K)‖∇a‖²` = the target. **Verdict: `K>0` fails on three independent counts; every positive route is H1.** **(a) Not weaker.** On flat `A` (`Ric=0`) weighted `CD(K,∞)` is *equivalent* to `Hess V_eff⪰K` (Bakry–Émery equality), so synthetic = pointwise = `G″⪰K` = H1 — "weaker than Hess" cannot exist in a flat chart; the sharp constant is the Savvidy bottom `K_*=infl(Hess S_W)+N/12=−g|B|+N/12<0`. **(b) Quotient gives no upgrade.** `A/G` is a compact-group quotient of a `CD(K_*,∞)` space; quotient stability *preserves* `K`, never raises it. The one apparent positivity — the O'Neill `A`-tensor *is* the field strength, `|A|²∝|F|²`, adding `+c|F|²≥0` to the base sectional — is a *sectional*, not the weighted-quotient Ricci, positive only where `F≠0`; at `F=0` (full-measure configurations) it vanishes and the Savvidy negative returns, while its **average** `E[c|F|²]` is the **gluon condensate `⟨¼trF²⟩`** — an H1 name. **(c) Uniformity.** `A/G` is compact (`D<∞`), so `CD(K_*,∞)`+Buser gives `ρ(K_*,D)>0`, but `ρ→0` as `D∼L²→∞`: a finite-`V` kinematic mass, not the gap. **Free massless is `RCD(0,∞)` (`K=0`, zero mode gapless, `ρ=0`) — acid PASS. #1,#4 PASS; #2,#3,#5 FAIL. NOT 5/5.**

## §1 The forge and the exact dictionary
`RCD(K,∞)` ⇒ heat-semigroup variance decay `d/dt Var(f_t)≤−2K Var(f_t)`, i.e. `Var(f_t)≤e^{−2Kt}Var(f)` (Savaré 2014; Gigli). For the coarse amplitude this is a *global* Poincaré `Var(a)≤K⁻¹ℰ(a)` at `K` — the correctly-signed object (upper bound on `Var` ⇒ lower bound on `G″`). Since `G″(k→0)=1/Var(a)` (D3, D2), `K>0` closes `m_eff²≥K`. The entire question collapses to **is `K>0` true, uniformly?**

## §2 On flat `A` the synthetic bound is the pointwise Hess — no weakening
`A≅ℝ^{(N²−1)V}` is flat, `Ric=0`; weighted measure `e^{−V_eff}dX`, `V_eff=S_W+Φ`. Bakry–Émery on a flat space: `CD(K,∞)` ⟺ `Hess V_eff⪰K` pointwise (Bochner equality). The Savvidy cone makes `Hess S_W` indefinite (`−D²+2adF_B`, negative transverse-magnetic mode, 7.1% post-quotient — real physics, W4); `Φ″⪰N/12` does not compensate (B3 §3, anti-localized). So `K_*=infl(Hess S_W)+N/12≈−g|B|+N/12<0`. There is no gap between synthetic and pointwise to exploit; and it is witnessed at a *smooth, μ-typical* point (the negative mode), so singularity-robustness — RCD's only edge over Hess — does not apply: the entropy is genuinely not `K`-convex for `K>0` along `W₂`-geodesics pushed in the negative direction.

## §3 Quotient: O'Neill positivity is the field strength, and its average is H1
`A/G` is singular (Gribov), where synthetic theory legitimately beats a nonexistent smooth Hess. But the quotient *preserves* `CD(K_*,∞)` (compact-group quotient never raises `K`), and the submersion `A→A/G` has `A_XY=−Ω(X,Y)` = field strength, so `sec_{A/G}=sec_A+¾|A|²=c|F|²` — the *one* positive term. It is **zero at `F=0`** (full-measure configs at strong coupling), where the Savvidy negative returns; global `K>0` needs it positive everywhere — false. Its mean `E[c|F|²]` *is* positive — the **gluon condensate `⟨¼trF²⟩`**, i.e. H1. So O'Neill either fails uniformity (pointwise) or restates the target (averaged): the frustration is exact.

## §4 The compactness loophole and its uniformity death
`A/G` compact ⇒ `D<∞`; Buser/Lott–Villani give `CD(K_*,∞)⇒PI(ρ)`, `ρ=ρ(K_*,D)>0`. But `ρ∼c(K_*)/D²→0` as `D∼L²→∞`: a finite-`V` kinematic mass — the lattice-`N/12` artifact (B3 §5), present at `β=0`, absent in the continuum. A lattice-spacing mass, not `m_eff²(k→0)≥λ` uniform in `(L,a)`.

## §5 Free-field acid test — clean PASS
Free massless `k→0`: zero mode flat, `p∝`Lebesgue over the (unbounded, spectral-dimension-reduced) field ⇒ Gaussian ⇒ `RCD(0,∞)`, `K=0`, `f₁²=0`, `G″=0`, `Var(a)=∞`, `ρ=0`. The synthetic bound correctly gives **no** gap and *would* be vacuous if finite — it is `K=0` exactly. Non-SUSY ✓. Pure YM's Savvidy `K_*<0` is the discriminator.

## Five-point scorecard
1. **Nonperturbative — PASS.** Global synthetic condition, all `β`, no `g`/`1/N`.
2. **Sign-controlling LOWER bound `m_eff²(k→0)≥λ` — FAIL.** `K>0` false: flat `A` ⟺ `Hess⪰K`, `K_*=infl(Hess S_W)+N/12<0`; witness = Savvidy negative mode (smooth, typical).
3. **Non-circular — FAIL.** On flat `A`, synthetic = `G″⪰K` = H1; on `A/G`, quotient preserves `K` (no upgrade); O'Neill positivity's average = gluon condensate `⟨¼trF²⟩` = H1; the `RCD⇒PI` constant *is* the coarse gap.
4. **Non-SUSY + fails free field — PASS.** `RCD(0,∞)`, `K=0`, zero mode gapless, `ρ=0`; exact.
5. **Uniform `(L,a)` — FAIL.** Only positive `ρ` is compactness (Buser), `ρ(K_*,D)∼1/D²→0` as `L→∞`; lattice artifact.

## Exact failing estimate
```
target:  G″(k→0) ≥ λ  ⟺  Var(a) ≤ 1/λ                [D3, H1]
RCD:     Var(f_t) ≤ e^{−2Kt}Var(f)  ⇒ PI λ₁ ≥ K        [Savaré 2014]
flat A:  CD(K,∞)  ⟺  Hess(S_W+Φ) ⪰ K  (Bakry–Émery equality)
         K_* = infl(Hess S_W)+N/12 = −g|B|+N/12 < 0    [Savvidy negative mode, 7.1%]
quotient: K preserved under compact-G quotient (never raised);
          O'Neill: sec_{A/G} = c|F|² ≥ 0  pointwise, =0 at F=0;
          E[c|F|²] = c⟨¼trF²⟩ = gluon condensate = H1   [circular]
compact: D∼L² , Buser ρ(K_*,D) ∼ c(K_*)/D² → 0 (L→∞)  [lattice artifact]
free:    RCD(0,∞), K=0, f₁²=0, Var(a)=∞, ρ=0           [acid PASS]
```

**Sharpest line.** Synthetic `RCD(K,∞)` is the most *geometric* of the lower-Ricci languages, but on the flat lattice chart it is **provably identical** to the pointwise Hess the Savvidy cone already falsifies (`K_*=−g|B|+N/12<0`) — there is no "weaker bound that survives." The gauge quotient's genuine new ingredient, the O'Neill `|F|²` curvature, is positive only where `F≠0`, vanishes on the `F=0` bulk, and averages to the **gluon condensate `⟨¼trF²⟩`** — the target wearing another name; the compactness Poincaré is a finite-`V` artifact vanishing as `L→∞`. Free massless is `RCD(0,∞)` with `ρ=0`, correctly failed. **The precise missing step: a *pointwise* (not averaged, not compactness) synthetic lower bound `K>0` on `(A/G, g, μ)` uniform in `(L,a)` — equivalently the positive-definiteness of the weighted Ricci on the quotient where `F=0`, which is verbatim H1. Passes 1,4; fails 2,3,5. NOT 5/5.**

*Refs.* (external) Lott–Villani, *Ann. of Math.* **169** (2009) 903; Sturm, *Acta Math.* **196** (2006) 65, 133; Ambrosio–Gigli–Savaré, *Duke Math. J.* **163** (2014) 1405; Savaré, *Discrete Contin. Dyn. Syst.* **34** (2014) 1641 (variance contraction / spectral gap); Ambrosio–Gigli–Savaré, *Ann. Probab.* **43** (2015) 339; Gigli, *J. Funct. Anal.* **267** (2014) 3450 (RCD spectrum); O'Neill, *Michigan Math. J.* **13** (1966) 459 (submersion / `A`-tensor); Bérard–Besson–Gallot, *Ann. Scuola Norm. Sup. Pisa* **9** (1982) 341; Buser, *Ann. Sci. ENS* **15** (1982) 213; Bakry–Émery 1985; Savvidy, *Phys. Lett. B* **71** (1977) 133. (internal) D3-SUSCEPTIBILITY-VARIATIONAL (`G″=1/Var(a)`); B3-FUNCTIONAL-BE (§3 anti-localized cone, §5 lattice-`N/12` caveat); PROOF-JACOBIAN-SAVES (`G″=Δ_dyn+W_geo+X`, Savvidy split); D2-COUPLING-DOEBLIN (`b′(0)=−G″`, von Renesse–Sturm `CD`); C8-SPECIES-NOGO (empty-referent / averaged-positivity pattern); OTHER-H1-POINCARE, PROOF-IR-H1 (H1⟺`m_eff²(k→0)≥λ`); CONTINUUM-VARIANCE-DEFECT §4 (anomaly ⇒ condensate existence, not heaviness — `⟨¼trF²⟩`=H1); PROOF-SAVVIDY-WGEO-DOMINATION (Savvidy indefinite, `ρ≈3–16`); k0-CONSTRAINT-EFFPOT (`V_eff″=G″`). No fabricated bound.
