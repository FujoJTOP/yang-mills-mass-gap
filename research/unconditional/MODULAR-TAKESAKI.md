# MODULAR-TAKESAKI — the vacuum modular operator is a mass-blind boost: a genuine free-field (entanglement-saturation) diagnostic that cannot cross to interacting YM, because its nonperturbative extension *is* the gap

> **判词 (verdict).** The species is forged cleanly and **passes the acid test in its only well-posed form** — the free *massless* vacuum has no modular/entanglement gap, so the diagnostic returns 0 there — but it **cannot lower-bound `m_eff²(k→0)` for YM** (fails #2, #3, #5). Two exact walls. **(analytic)** On the continuum wedge the modular Hamiltonian is, by Bisognano–Wichmann, the boost `K_W=U(Λ_W)`; its spectrum is all of ℝ and **independent of the mass**. The premise "BW geometric flow *encodes* the gap" is false: `spectrum(Δ_W)` is continuous and identical for `m=0` and `m>0`, so every functional of the wedge modular operator is **mass-blind** — 0 for free massless *and* free massive alike. **(structural)** The one non-vacuous incarnation, the finite-region **entanglement spectrum** (Li–Haldane), *does* see the gap in **Gaussian** systems (entropy saturation ⟺ `m>0`, confirmed below), but extending "saturation ⟺ gap" to interacting 4D YM is exactly the Li–Haldane-type conjecture, and controlling the interacting modular Hamiltonian **is** the gap. The Araki relative-entropy route `S(ω_glueball‖ω_vac)` is no better: UV/area-law-dominated (→ free by asymptotic freedom), and it needs a "glueball state" that is *defined* by the spectrum. **No closure; #1 genuine, #4 only vacuous.**

## 1. The target is a spectral-mass object; the modular operators are geometric
The narrowed target `E[G′]−Var(G)=G″≥λ` is the stiffness of the volume zero-mode `a₀=⟨v,X⟩`: `G″=V_eff″(0)=1/Var(a₀)` (k0-CONSTRAINT-EFFPOT; PROOF-JACOBIAN-SAVES). By the spectral theorem `Var(a₀)=∫|⟨Ω|a₀|λ⟩|² dμ_H(λ)` over the **physical Hamiltonian** spectrum, and `m_eff²(k→0)` is the infimum of the support of the `k=0` spectral measure (Källén–Lehmann threshold). This is an **IR, Hamiltonian-spectral** object. The modular data — `Δ_W`, `K_W=−logΔ_W`, the relative entropy — live on the **algebra of a region** and, for the vacuum of any Poincaré-covariant theory, are controlled by covariance + Reeh–Schlieder + wedge duality, *not* by the mass spectrum.

## 2. The acid test — passed, but it exposes the tool as free-only
For a free field the finite-region modular Hamiltonian is the entanglement Hamiltonian, and its entropy is a faithful Gaussian diagnostic. Half-chain entropy `S` of a 1D free-fermion chain (numpy, L up to 800), staggered mass `m`:
```
 m=0.0 :  S = (1/6)ln L + 0.287  (log-divergent — massless, NO gap)
 m=0.1 :  S = 0.7292 (L-independent)   m=0.4 : S = 0.4887 (L-independent)
```
So the diagnostic gives **no saturation at `m=0`** (`f₁²=0`, acid test → 0, PASS #4) and **saturation for `m>0`**. It is a correct free-field gap detector — *precisely the Gaussian regime*. It is not a YM tool.

## 3. The analytic wall — Bisognano–Wichmann is mass-blind
`K_W=U(Λ_W)` for the wedge, in *any* theory satisfying the Wightman axioms. The boost generator's spectrum is ℝ (continuous) and does not depend on `m`; the free massless and free massive vacua have the **same** `spectrum(Δ_W)`. Hence a "modular spectral gap" of `Δ_W` is identically absent, and any `λ` it produced would be `0` for all masses — it cannot give `λ>0` for YM (#2, #5 fail). Half-sided modular inclusions (Borchers) recover the *physical translation*, and hence a gap, **only when the modular flow is geometric** — free, conformal, or dilation-covariant — i.e. the regime complementary to interacting YM; the Borchers identification is open/false otherwise.

## 4. The structural wall — relative entropy and interacting ES are the gap relocated
(a) **UV-dominated:** region entanglement entropy is area-law-divergent; asymptotic freedom makes YM's UV free, so UV modular data cannot see the `k→0` threshold. (b) **Interacting ES is conjectural:** Li–Haldane "ES gap ⟺ gap" is proven only for free/Gaussian systems; for 4D YM the interacting modular Hamiltonian is intractable. (c) **Definitional circularity:** `S(ω_glueball‖ω_vac)` presupposes "the glueball," *defined* as the first `0⁺⁺` excitation — the gap. (d) **Bridge = RG:** mapping modular data to zero-mode stiffness needs the RG flow CLAUDE.md §1 says does not exist.

## 5. Exact missing step
```
∃λ>0:  F[Δ_W, S_rel] = 0 (free massless)  ∧  F ≥ λ (YM)  ⟺  G″(k→0) ≥ λ.
```
Wedge `F` is mass-blind ⇒ `F≡0`. Finite-region `F` needs "entanglement saturation ⟹ Hamiltonian spectral gap" for interacting YM — open (Li–Haldane), i.e. **F = the gap**. The species is a RENAME of the (CP)/H1 wall in modular language, differing only in *why* it fails: not torsion, not a false sign (C8), but a **mass-blind operator matched to an IR spectral target**.

## Five-point verdict (honest)
1. **Nonperturbative — PASS.** Modular theory / relative entropy is exact; no loop expansion.
2. **Sign lower bound `m_eff²≥λ` — FAIL.** Wedge `Δ_W` is a mass-blind boost; no positive `λ`.
3. **Non-circular — FAIL.** Only non-vacuous form (interacting ES saturation ⟺ gap) *is* the gap; relative entropy is controlled by vacuum modular data whose IR is the gap.
4. **Non-SUSY + fails free field — PASS, vacuously.** Free massless → 0 confirmed (log-divergent `S`, `f₁²=0`), but the wedge form also gives 0 for free *massive*; non-vacuous form is free-only.
5. **Uniform lower bound — FAIL.** No uniform `λ`; `0` for all masses.

**Sharpest line.** The species names a real free-field object — entanglement-spectrum saturation — and passes the acid test by *being* that free-field object. Two things kill it: the continuum wedge modular Hamiltonian is the boost, spectrally independent of `m` (so "modular spectral gap" is the wrong operator for a Hamiltonian gap), and the only bridge from modular/entanglement data to the zero-mode stiffness `G″(k→0)` is the interacting RG flow — the gap itself. Modular theory describes the *geometry* of the vacuum, not the *stiffness* of its softest mode.

*Refs.* Bisognano–Wichmann, J. Math. Phys. **17** (1976) 303; Araki, Publ. RIMS **11** (1976) 809 (relative entropy); Borchers, CMP **143** (1992) 327 (half-sided modular); Wiesbrock, CMP **157** (1993) 83; Li–Haldane, PRL **101** (2008) 010504; Casini–Huerta; Blanco–Casini–Hung–Wood (relative-entropy/Bekenstein bound); Reeh–Schlieder. Internal: k0-CONSTRAINT-EFFPOT (`V_eff″=G″=1/Var(a₀)`); INVENT-VARIANCE-INEQ (`G″=S−D`); PROOF-JACOBIAN-SAVES (`G″=Δ_dyn+W_geo+X`); PROOF-IR-H1 (H1 = gap); PROOF-GAUSS-INVARIANT-PF (free Maxwell: PF vacuum, yet `f₁²=0`); OTHER-H1-POINCARE; C8-SPECIES-NOGO.
