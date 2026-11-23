# PIROGOV-SINAI — the screen-cluster expansion is exact, passes the acid test, and still does not close: the Debye mass *is* the k→0 susceptibility = `1/Var(a)` = H1, its convergent radius is `P3`'s `β<β₁`, and its only extension forfeits the sign

> **判词 (verdict).** The forged species — a Pirogov–Sinai contour / Mayer-resummed Debye-screening cluster expansion giving exponential decay of the 0⁺⁺ correlator at all β below the transition — does **not** close. **(i)** Where it is a genuine absolutely-convergent bound (`β<β₁=O(1/N)`), the screening mass *equals* the KP/Osterwalder–Seiler result: **P3 in new clothes**. **(ii)** Extending past `β₁` needs a Borel/analytic continuation, **not sign-controlling**, carrying the gluon-condensate ambiguity. **(iii)** Throughout, the screening mass `κ_D²=lim_{k→0}G(k)` is the k→0 susceptibility, which by the Ward/Stein identity `G″·Var(a)=1` is exactly `1/Var(a)` = **H1** — a restatement, not an estimate. Free-field acid test **PASSES** (genuine: `κ_D∝` interaction, `Var(a)=∞`) but is necessary, not sufficient. Five-point: **1 PASS (regime-limited), 2 FAIL, 3 FAIL, 4 PASS, 5 FAIL ⇒ 2/5**; analyticity wall untouched.

## §1 The forge — what a lower bound actually requires
The 0⁺⁺ channel is sourced by the plaquette density `ρ=(1/N)Re tr F²`, `G(x)=⟨ρ(x)ρ(0)⟩_c=Σ_{n≥1}|c_n|²e^{−Δ_n t}`. Debye–Hückel chain resummation gives `G(k)=κ_D²/(k²+κ_D²)·(bare)`, `κ_D²=lim_{k→0}G(k)=density susceptibility`. A **lower** bound `m≥λ` needs **two** facts (as in P3): **(a)** an exponential upper bound `|G(t)|≤Ce^{−λt}`, **and (b)** a residue lower bound `f₁²≥r₀>0`; the spectral step `r₀e^{−Δ₁t}≤G(t)≤Ce^{−λt}` then forces `Δ₁≥λ`. A screen-cluster species is honest only if it supplies *both* from first principles, uniformly in `(L,a)`. Neither is free past `β₁`.

## §2 Two independent deaths
**(A) Channel/phase mismatch.** T=0 confined YM is **not** a Debye plasma: its long-distance physics is a confining *string* (area law), and the lightest 0⁺⁺ is a closed flux tube / glueball, not a screened plasmon. The "screening mass" YM actually has is either the dual-superconductor Meissner mass of a 2-form (spin-1 scale `√σ`, **cross-channel**, H1-T10-dual) or a finite-T `m_D~gT` (thermal, ≠ T=0 scalar). The forge resums the wrong object in the wrong phase.
**(B) Circularity / analyticity dichotomy.** `κ_D²=lim_{k→0}G(k)` is the very correlator to bound; by `G″·Var(a)=1`, `κ_D²≥λ² ⟺ Var(a)≤1/λ²` = H1 — an *identity*, not an estimate. Making it a bound needs a small bare activity (`u≍β`) ⇒ `β<β₁`, reproducing P3; crossing needs analytic continuation, but the absolute-convergence radius is `β₁` and, though the mass *is* β-analytic across the smooth crossover (PROOF-NO-MELT-PHASE), the absolute bound cannot be continued (ATOM-strong-coupling-extend) and a Borel continuation is not sign-controlling (ambiguity = gluon condensate). Sharper, PS-specific: contours need a **stable reference** with small activities, but the Savvidy one-loop Hessian `−D²+2adF` is indefinite (transverse magnetic `−g|B|<0`, W4 7.1% negative) — **no perturbative stable vacuum** exists around which to define contour activities in the intermediate regime.

## §3 Free-field acid test — the one genuine pass
A correctly built screen expansion has `κ_D²∝interaction`: the free massless Gaussian has a field-independent determinant (no charge density to resum), so `κ_D=0`, `G(x)~` power-law, `Var(a)=∞`, `f₁²=0`, and the resummation correctly yields **no mass** — it does not fabricate a gap. At `β<β₁`, `u≍β>0 ⇒ κ_D>0 ⇒` exponential. The tool genuinely distinguishes screening YM from the free Coulomb field; this is its real merit and why it is **not vacuous** (unlike B4 route (b)). But the acid test is the `β=∞` boundary; it says nothing about the wall at `β≈2.5`, where the activity is `O(1)` and the convergent resummation *fails to converge*.

## §4 Five-point
1. **Nonperturbative — PASS, regime-limited.** Contour/PS and Mayer resummations are exact (no loop truncation) at `β<β₁`; the only Borel extension is perturbative-asymptotic.
2. **Sign-controlling LOWER bound `m_eff²(k→0)≥λ` — FAIL.** Inside the radius it holds (=P3); past it the residue bound `f₁²≥r₀` is exactly the β-nonanalytic-jump no-go (PROOF-F1SQ), and the variance form collides with equipartition `Var(G|0)≈2 ⇒ G″≈E[H]−Var≈0`, never `≥λ` (B4).
3. **Non-circular — FAIL.** `κ_D² ≡ k→0 susceptibility ≡ 1/Var(a) ≡ H1` (Ward); the resummation re-expresses, does not bound.
4. **Non-SUSY + fails free field — PASS (genuine, §3).**
5. **Uniform `≥λ` in `(L,a)`, all β — FAIL.** Screening radius `β₁`; continuation non-uniform (Borel ambiguity); `a→0` open.

## §5 Exact failing estimate
```
target:  G″(0)≥λ ⟺ Var(a)≤1/λ            (Ward G″·Var(a)=1)
Debye:   κ_D² = lim_{k→0}⟨ρρ⟩ = 1/Var(a) = H1   (identity, not bound)
radius:  u(β)K<1 ⇒ β<β₁=O(1/N)             (=P3; ATOM/screening wall)
residue: f₁²>0 across [1/48,2.5] = β-nonanalytic-jump no-go  (F1SQ)
equipart: Var(G|0)≈2 ⇒ G″≈E[H]−Var≈0, never ≥λ   (B4)
free:    κ_D=0, Var(a)=∞, f₁²=0 ⇒ correctly no mass (acid PASS)
Savvidy: Hess=−D²+2adF indefinite (−g|B|<0) ⇒ no stable PS reference
```

**Sharpest line.** The screen-cluster expansion is exact and passes the acid test — which is exactly why it cannot help: its screening mass is the zero-momentum susceptibility it is meant to bound (`κ_D²=1/Var(a)=H1`, an identity), its convergent radius is the strong-coupling radius `β₁` it must cross (P3), and its only continuation is a non-sign-controlling Borel sum whose ambiguity is the gluon condensate. YM at T=0 is a confining string, not a Debye plasma, so the screening it "resums" is either the cross-channel Meissner `√σ` (H1-T10-dual) or a thermal `m_D`. The species contributes a genuine free-field discriminator and **zero reach past the analyticity wall: 2/5, wall untouched.**

*Refs.* Pirogov–Sinai 1975/76; Zahradník (PS reviews); Kotecký–Preiss, CMP **103** (1986) 19; Brydges–Federbush (screening / dipole-gas resummation); Brydges–Giffen–Durhuus–Fröhlich, NPB **275** (1986) 459; Osterwalder–Seiler, Ann. Phys. **110** (1978) 440; Shen–Zhu–Zhu, CMP **400** (2023) 805, arXiv:2204.12737; Cao–Nissim–Sheffield, arXiv:2505.16585 / 2509.04688; Polyakov, NPB **120** (1977) 429 (3D Debye=gap, wrong `d`). Internal: PROOF-P3-UNCONDITIONAL; PROOF-P1P2; PROOF-F1SQ; ATOM-strong-coupling-extend; B4-VARIANCE-MAG; B5-CORRELATION-CONTROL; C7-SPECIES-CONSTRUCT; C8-SPECIES-NOGO; H1-T10-dual; PROOF-CONFINEMENT-RANDOM-SURFACE; PROOF-NO-MELT-PHASE; G1-UNIFORM-ZA; PROOF-JACOBIAN-SAVES; w7 measures (`VarG≈2`, `W_geo≈0.32`).
