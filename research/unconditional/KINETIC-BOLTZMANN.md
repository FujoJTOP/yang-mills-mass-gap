# KINETIC-BOLTZMANN — the collision gap is the fourth lens on `G″`; the physical Boltzmann frequency is a sign-positive sum of squares that bounds thermalization, *not* the mass gap, and the Savvidy tachyon blocks the strict floor

> **2026-11-23 · NEW SPECIES (kinetic / Boltzmann / relaxation).** 判词 (verdict). Two collision operators exist for Yang–Mills; they fail oppositely. **(A) The exact linearized generator** (Fokker–Planck/Lindblad) is a genuine non-static species, but its gap on the coarse mode is the drift curvature `−β′(0)=G″` by construction — the *same* spectral number met as the `W₂`-rate (D2), the Bakry–Émery/LSI constant (E1), and the Krein–Rutman gap (F2). Four lenses, one number; von Renesse–Sturm is sharp, so `gap≥λ` *is* `G″≥λ` is H1. **(B) The physical two-body Boltzmann operator** gives a collision frequency `ν ≃ Σ g⁴|M|² ≥ 0` — a *sum of positive terms*, the one structural advantage over every Hessian route — but `ν` bounds **thermalization of hard modes**, not the `k→0` pole; the long-wavelength sector of `C_B` is hydrodynamic (`∼Dk²→0`), and at `T→0` the only kinetic scale `gT→0`. Neither yields `G″≥λ`. **Free-field acid test PASSES cleanly** (`g=0` ⇒ `C≡0`; `k=0` flat ⇒ rate 0). **#1,4 PASS; #2,3,5 FAIL. 2/5 — no closure.**

## §1 Two collision operators
Equilibrium `μ∝e^{−S_W}dU` on `𝒜=G^E`; coarse amplitude `a=⟨v,X⟩` (`k→0` block mode). **(A)** Linearize the exact reversible semigroup (classical FP; the classical limit of a Lindblad generator — `LINDBLAD-DECOHERENCE.md` is absent, but the structure survives): a self-adjoint `C_A` with detailed balance and an H-theorem. **(B)** Two-body collisions among gauge excitations, `C_B f=∫(f'f'_1−ff_1)|M|²` — effective kinetic theory (Arnold–Moore–Yaffe). Both non-static; both fail.

## §2 Exact generator — its gap *is* `G″` (one object, four lenses)
For `C_A` on `a` the relaxation rate is `−β′(0)=G″₀=S−D` (E1 §1, D2 §1). Not a new bound: `gap_{C_A}` = the `W₂`-rate (D2, von Renesse–Sturm) = the LSI/Bakry–Émery constant of `V_eff` (E1) = the Krein–Rutman top eigenvalue (F2). `gap = W₂-rate = LSI = operator-gap = G″ = S−D`. von Renesse–Sturm is an *equivalence*, so nothing kinetic beats it: `gap≥λ ⟺ G″≥λ` ⟺ H1. **Circular (#3).**

## §3 H-theorem — existence of relaxation, not a rate
Entropy production `dH/dt=−E[(∇log ρ/μ)²]≤0` (Fisher information) proves the semigroup relaxes to `μ` nonperturbatively — but its magnitude **degenerates on the coarse/flat directions**: at `k→0` the score vanishes along the Savvidy-negative direction, so no *lower* rate there. The H-theorem is a second law (existence), not a Poincaré bound (rate); turning it into a rate is the LSI = §2 = H1.

## §4 Physical Boltzmann — a positive sum of squares bounding the wrong quantity
`C_B` is the honest "collision frequency": its gap `ν ≃ Σ g⁴|M|² ≥ 0` is manifestly *positive at any coupling* — no indefinite Hessian to fight. But it is the wrong positivity: `ν` is the rate at which *hard* quasi-particles thermalize, UV-dominated (`∝g⁴T`), whereas the mass gap is the *slowest* physical pole. `ν≫m` says nothing about the slow glueball; at `T→0` the only kinetic mass `gT→0`. Worse, the `k→0` sector of `C_B` is **hydrodynamic** (`ω∼Dk²→0`): the collision gap tends to `0` on the very mode the target asks about. `ν≥λ` is true and irrelevant — it does not imply `G″≥λ` (**#2**) — and effective kinetic theory is a `g≪1` quasiparticle expansion, so `#1` is only nominal for B.

## §5 Acid test (clean)
Free massless: no vertices ⇒ `C_B≡0` (no collisions, no relaxation); `k=0` coarse mode is a flat direction (`⟨v,∇S⟩≡0`) ⇒ `C_A`-gap `=0`. The tool returns `λ=0` — it does *not* manufacture a gap. Massive free: rate `m²`. Pure Haar `β=0`: `W_geo≥N/12` (Prékopa; PROOF-JACOBIAN-SAVES). Exact ⇒ non-vacuous (**#4 PASS**).

## Five-point scorecard
1. **Nonperturbative — PASS (A).** Exact FP/Lindblad generator, no loop expansion, all `β`; B is weak-coupling (nominal).
2. **Sign-controlling `G″≥λ` — FAIL.** A circular; B bounds `ν` — positive but `≠G″`, `→0` at `k→0`.
3. **Non-circular — FAIL.** `gap_{C_A}=G″` exactly; four lenses, one number.
4. **Non-SUSY + fails free field — PASS (cleanest).** `C≡0`, flat `k=0` ⇒ `λ=0`; massive `m²`; Haar `N/12`.
5. **Uniform `≥λ` — FAIL.** `ν` UV-dominated/non-uniform; no uniform floor on `G″`.

## Exact failing estimate
```
Goal:  gap(C) ≥ λ  (⇒ G″(k→0) ≥ λ).
(A) C = exact linearized generator (FP/Lindblad):
      gap_C = −β′(0) = G″(0) = S − D = Hess V_eff
            [von Renesse–Sturm / Bakry–Émery / Krein–Rutman — SHARP]
      ⇒ "gap ≥ λ" ⟺ G″ ≥ λ  (H1).                         [CIRCULAR]
(B) C = two-body Boltzmann collision operator:
      gap_C = ν = Σ g⁴|M|² ≥ 0     (H-theorem; sign-positive, any coupling)
      but ν bounds thermalization of HARD modes, not the k→0 pole;
      k→0 sector hydrodynamic (ν_hydro=0); ν ~ g⁴T → 0 at T→0.
      ⇒ ν ≥ λ  does NOT give  G″ ≥ λ.                     [WRONG QUANTITY]
free massless: C ≡ 0 (no vertices), k=0 flat ⇒ gap 0.     [ACID PASS]
```

**Sharpest line.** The kinetic lens supplies what nobody else does — an *exact nonperturbative generator* whose gap is a spectral number, and a *manifestly positive* bracket `Σ|M|²≥0` that sidesteps the Savvidy indefinite Hessian — and both miss the target for one reason. Route (A)'s gap is literally `G″=S−D`, the fourth face of H1 (`W₂`=LSI=KR=kinetic-gap), sharp and un-improvable; route (B)'s positivity is real but lives on *hard* UV modes, while the mass gap is the slowest IR pole and the `k→0` collision spectrum is hydrodynamic (`→0`). "A positive sum of squares replaces the indefinite Hessian" dies because that sum bounds the thermalization rate, not the gap: fast collisions coexist with a gapless long-wavelength sector. The strict floor at the Savvidy crossover is still H1; the kinetic species renames it and points its positivity at the wrong observable. **#1,4 pass; #2,3,5 (the Clay content) fail. No fabrication.**

*Refs.* Boltzmann, Sitzungsber. Akad. Wiss. Wien **66** (1872) 275; Cercignani, *The Boltzmann Equation and Its Applications*, AMS 67, Springer (1988); Cercignani–Illner–Pulvirenti, *The Mathematical Theory of Dilute Gases*, AMS 106, Springer (1994); Guo, *Indiana Univ. Math. J.* **53** (2004) 1081 (spectral gap of linearized operator); Arnold–Moore–Yaffe, *JHEP* **0011** (2000) 001 [hep-ph/0005074]; *JHEP* **0305** (2003) 051 [hep-ph/0302165]; Lindblad, *CMP* **48** (1976) 119; Gorini–Kossakowski–Sudarshan, *J. Math. Phys.* **17** (1976) 821; Wilson, *PRD* **10** (1974) 2445; Lüscher, *CMP* **54** (1977) 78; von Renesse–Sturm, *CPAM* **58** (2005) 923; Bakry–Émery, *LNM* **1123** (1985) 177. Internal: E1-FRG (exact FP generator); D2-COUPLING-DOEBLIN (`b′(0)=−G″`); F2-KREIN-RUTMAN-QUANT (operator gap = `G″`); C8-SPECIES-NOGO (no sign-control of cross terms); INVENT-VARIANCE-INEQ (`G″=S−D`); PROOF-COARSE-POINCARE (`PᵀHP≈−0.3`); PROOF-JACOBIAN-SAVES (`W_geo≥N/12`); W7-Gpp0 (`G″=+0.51`). No fabricated bound.
