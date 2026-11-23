# G2E-RIGOROUS-NU — on the rigorous KL measure the gap is the *support* problem; the sufficient target `inf supp ν>0` factorizes exactly as (residue `f₁²`)×(susceptibility `Var a`), and P3 supplies only the first factor, only at β<1/48

> **判词 (verdict).** Grant the OpenAI/Lean-2026 premise, hence a genuine Wightman KL measure `ν` with `ν≥0`, `supp ν⊂[0,∞)` (G1-CONTINUUM-OS). The task's retarget is **correct**: the Clay datum is `Δ>0 ⟺ inf supp ν=Δ²>0`, and the campaign's `H1 = Var(a)≤1/λ` is a **necessary-but-not-sufficient** proxy (REDUCTION-AUDIT §(b); G2H §2). On the rigorous `ν` the two are tied by one elementary spectral inequality — `χ := Var(a) = ∫ s⁻¹dν(s) ≥ f₁²/m₁² ⇒ m₁² ≥ f₁²/χ` — so the **sufficient** data is the *pair* `(f₁², Var a)`, never `Var a` alone. Running the proven pieces on `ν` confirms this: P3 gives the residue `f₁²≥r₀` only at β<1/48; `W_geo≥0` bounds only the *sign* of the geometric term; WSTC bounds the wrong (high-field) region; no-dissolution is conditional on the gap. **They do not assemble to `inf supp ν>0` at β≈2.5**, and `H1` alone can never (the tiny-residue counterexample is rigorous on the KL measure). Not 5/5; Clay 0%.

## §1 · The correct target, and why retargeting also removes the Ward non-exactness
Write `m₁² := inf supp ν`, `f₁² := ν({m₁²}) = |⟨0|trF²|G⟩|²`, `M₀ = ∫dν`, `χ = ∫s⁻¹dν`. Because the atom at `m₁²` contributes `f₁²/m₁²` and every other mass `s≥m₁²` contributes `1/s≤1/m₁²`, one has the **exact** measure identity
```
   χ = f₁²/m₁² + ∫_{[m₂²,∞)} s⁻¹ dν  ≥ f₁²/m₁²   ⟹   m₁² ≥ f₁²/χ .      (Bridge)
```
Hence `inf supp ν>0` is implied by `(f₁²≥r₀) ∧ (χ≤1/λ)`, giving `m₁² ≥ r₀λ>0` and a channel gap at rate `√(r₀λ)`. This **sidesteps** the task's `G″(0)=1/Var(a)`: that holds only for a Gaussian marginal (`G″·Var(a)=1` iff Gaussian; 2–5% off measured — ABC-NEW-INEQ §1, REDUCTION-AUDIT §(C)). The Bridge is exact for *any* positive measure — no Gaussianity — so retargeting corrects H1's insufficiency **and** discards the only inexact step. By KL (G2H P1) `inf supp ν=Δ²` is *equivalent* to uniform exponential clustering: the true, not the lossy, condition.

## §2 · The pieces, run on the rigorous `ν`
- **P3 residue** (PROOF-P3-UNCONDITIONAL): `f₁²≥r₀>0` at β<1/48, V,L-uniform — **factor 1**, genuine but regime-limited.
- **`W_geo≥0`** (Prékopa, PROOF-JACOBIAN-SAVES): `G″=Δ_dyn+W_geo+X` fixes the **sign** of the geometric term on every mode incl. `k→0`, but `G″≥λ` still needs an upper bound on `Var(G|0)`=`χ` — so it supplies **neither** factor.
- **WSTC** (ABC-NEW-INEQ §4): high-field cone variance `≤e^{−cβ(κ²−κ_typ²)N_□}→0`, but sign-carriers sit at `κ_typ` — it bounds where `ν` has **no spectral mass**; silent on the bottleneck.
- **No-dissolution** (kinematic): protects an atom below `2m_g` **given** a uniform gap + `0⁺⁺`-lightest ordering (both the species): it preserves, not produces.

## §3 · Does the well-posed problem assemble to a gap? No — and here is exactly why
At **β<1/48** the assembly *does* close: SZZ/P3 give `f₁²≥r₀` **and** LSI⇒finite `χ`, so the Bridge yields `m₁²≥r₀λ>0` — but this merely **reproduces** the known strong-coupling theorem (confirming the Bridge is correct, adding nothing). At **β≈2.5** it fails on **both** factors: `f₁²` is unproven, and transport `f₁²(μ)=Z_O²f₁²(μ₀)` is exact but **presumes `(CP)`=single phase=`I1∧I2`=the gap** (circular); `H1` is open at every β. The rigorous object does not weaken either — *same* measure, now well-defined. **Gauge-survival caveat:** a scaling limit that loses the SU(N) group is plausibly free/massless (`inf supp ν=0`), so "`ν` is available" is itself conditional on the property neighboring the gap.

## §4 · Acid test (free massless field)
`f₁²=0`, `χ=∫s⁻¹dν=∞` ⇒ Bridge gives `m₁²≥0` (vacuous, *correct*), `inf supp ν=0`; `W_geo`, WSTC, no-dissolution all return nothing. **No spurious positivity — PASS**, on every piece.

## §5 · Honest 5-point
1. **Nonperturbative — PASS.** Exact measure identity, all β, no `1/N`/`g`.
2. **Sign-controlling — PARTIAL.** The Bridge is a right-sign *lower* bound on `m₁²` **given** `f₁²`; but there is no β≈2.5 `f₁²`, so no bound there.
3. **Non-circular — FAIL at closing.** Residue transport across `[1/48,2.5]` requires `(CP)` = the gap.
4. **Non-SUSY + fails free field — PASS (clean).**
5. **Uniform `inf supp ν>0` at β≈2.5 — FAIL.** Needs **both** factors; neither is available at 2.5.

**Exact missing step.** A uniform-in-`(β∈[1/48,2.5],a,V)` lower bound `f₁²≥r₀>0` on the `0⁺⁺` KL **residue** (or a non-circular proof of `(CP)` enabling its transport) — the **residue factor** — *together with* `H1`. On the rigorous `ν` the gap is provably a **product**, and no variance bound, `W_geo≥0`, WSTC, or no-dissolution supplies the residue factor. That residue datum is the one thing the rigorous object makes *sharply nameable* but does not provide.

**Sharpest line.** The rigorous `ν` does not create a mechanism; it makes `inf supp ν>0` well-posed and acid-clean, and reveals that the campaign's "one species `H1`" is a *necessary* slogan: the *sufficient* object is the ordered pair `(residue, susceptibility)`, glued by the exact Bridge `m₁²≥f₁²/χ`. P3 furnishes the residue only at strong coupling; everything else is sign, wrong-region, or conditional. Clay stays a conjunction of a residue bound at β≈2.5 and an H1 susceptibility bound — **one tool form, two factors, zero closed at 2.5; 0%.**

*Refs (internal).* G1-CONTINUUM-OS (object≠gap; `Δ>0⟺inf supp ν`); REDUCTION-AUDIT §(a)–(c) (H1 necessary-not-sufficient; susceptibility≠mass; tiny-residue counterexample); G2H (susceptibility trap; gap⇔uniform clustering; (MS)); G2G-TRANSFER-SPECTRAL (rate=`G″`=H1, isolation=H*); G2-ASSEMBLE-100 §2,§7; GAP-G2-edge (residue–level-gap; `I1⟺I2`); PROOF-F1SQ (gapped free `:φ²:` has `f₁²=0`); PROOF-JACOBIAN-SAVES; ABC-NEW-INEQ §1,§4; PROOF-P3-UNCONDITIONAL; PROOF-NO-DISSOLUTION-KINEMATIC; A4-SPECTRAL-BOOTSTRAP (moments blind to `inf supp`). *External.* Osterwalder–Schrader, CMP **31** (1973) 83, **42** (1975) 281; Osterwalder–Seiler, Ann. Phys. **110** (1978) 440; Källén 1952; Lehmann 1954; Shen–Zhu–Zhu, CMP **400** (2023) 805, arXiv:2204.12737.
