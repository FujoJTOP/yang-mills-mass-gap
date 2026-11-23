# HYDRODYNAMIC-GAP — the long-wavelength EFT of the `0⁺⁺` sector is the block-spin diffusion itself; the "relaxation gap" `Γ₀` is verbatim `G″(k→0)=m_g`, so "absence of a diffusive mode" is the gap restated — forced modes are mass-blind, and the one discriminating mode is the trace anomaly (a scale-detector, not a gap-detector)

> **2026-11-23 · NEW SPECIES — long-wavelength hydrodynamic / effective field theory.** 判词 (verdict). Forge the hydrodynamic EFT of the `0⁺⁺` sector (Mori–Zwanzig slow-mode reduction on the gauge-invariant scalar `Tr F² = θ^μ_μ/(β/2g)`) and ask whether the **absence of a gapless diffusive mode / the presence of a relaxation gap `Γ₀`** bounds `m_eff`. **Verdict: this is a language, not a lever.** The `0⁺⁺` collective mode *is* the block-spin coordinate `a` (D2), whose Mori–Zwanzig relaxation rate at `k→0` is **identically `G″(k→0)=Γ₀`**; "relaxational (gap `Γ₀>0`) vs diffusive (rate `∝k²→0`)" *is* the gapped/gapless dichotomy, so demanding `Γ₀≥λ` **is** H1. The conserved-sector diffusive modes (sound, shear) are **conservation-mandated and mass-blind** — present in free massless, free massive, and gapped YM alike, so their presence cannot signal gaplessness (C8/AQFT empty-referent). The free field *is* correctly failed (`T^μ_μ=0` conformal ⇒ scalar conserved ⇒ gapless diffusive ⇒ `Γ₀=0`), but the discriminator is the **trace anomaly** (`β(g)≠0`, a *scale*), which certifies "a mass exists," never "`Γ₀≥λ`." **#1,#4 PASS; #2,#3,#5 FAIL = the whole gap. No `λ`. NOT 5/5.**

## §1 The forge and the exact dictionary
The `0⁺⁺` sector is carried by gauge-invariant scalars; the lightest is `Tr F²`, the trace-anomaly density `θ^μ_μ=(β/2g)Tr F²`. Mori–Zwanzig on the slow amplitude `a=⟨v,X⟩` (`k→0` block mode, D2) gives the hydrodynamic effective theory
```
∂_t A = iΩA − ∫_0^t K(t−t′)A(t′) ,   A ∈ 0⁺⁺  scalar.
```
Poles in Laplace space classify the collective mode:
- **conserved density** (`∂_t n+∇·j=0`, the energy/shear sector of `T^{μν}`): Fick ⇒ `ω=−iDk²` — **gapless, diffusive**;
- **non-conserved scalar** (the glueball `Tr F²`): no `iΩ` coupling to a current, finite memory `K(0)` ⇒ `ω=−iΓ₀`, with `Γ₀=(1/χ)∫_0^∞dt⟨Ė Ė⟩_{disc}>0` — **relaxation gap**.

The **hydrodynamic dictionary**, read against D2 (`b(a)=G″a`, `b′(0)=−G″`): the coarse mode's `k→0` relaxation rate *is* the curvature,
```
Γ₀ = −b′(0) = G″(k→0) = S−D = 1/Var(a) .          (identity, not estimate)
```
Hence `G″(k→0)≥λ` ⟺ the `0⁺⁺` collective mode has a relaxation gap `Γ₀≥λ` ⟺ `Var(a)≤1/λ` = **H1**, verbatim.

## §2 Why "absence of a diffusive mode" is the gap restated
A mode is **diffusive** iff its `k→0` rate `→0` (`∝k²`), i.e. gapless; **relaxational** iff the rate `→Γ₀>0`, i.e. gapped. So "the coarse mode has a relaxation gap `Γ₀`" and "`G″(k→0)=Γ₀>0`" say the **same thing** — the hydrodynamic EFT *classifies* the mode by the dichotomy it cannot *decide*. To compute `Γ₀` you need `∫_0^∞dt⟨ĖĖ⟩`, i.e. the `ω→0` spectral density of `Tr F²`, i.e. the pole position `m_g` — the gap. **Circular at line one.**

## §3 The conserved modes are mass-blind (the acid-test core)
The genuinely *hydrodynamic* (diffusive/sound) modes live in the conserved `T^{μν}` sector, protected by `∂_μT^{μν}=0`. They exist in **every** relativistic QFT — free massless, free **massive** (`Γ` the shear rate), and gapped YM — with transport coefficients `(D,η,c_s)` fixed by thermodynamics, not by the gap. Therefore their **presence cannot signal gaplessness** and their **absence** in a non-conserved scalar is automatic. This is the C8 empty-referent / AQFT mass-blindness pattern in kinetic dress: a mass-insensitive class, no inequality in it separating gapped from gapless.

## §4 Free-field acid test + the one real discriminator
**Free massless:** `T^μ_μ=0` (conformal, `β=0`), so `Tr F²` is (a piece of) the *conserved* stress tensor ⇒ the gapless diffusive (sound/shear) mode is **present**, `Γ₀=0`, `G″=0`. The tool **correctly gives no gap** — acid test **PASSED** (not the AQFT-style refutation, where the free field satisfies the hypothesis yet is gapless). **Gapped YM:** `θ^μ_μ=(β/2g)Tr F²≠0` ⇒ the scalar is **non-conserved** ⇒ relaxation, no diffusive pole. So free vs YM *are* separated — **by conformality (trace anomaly), not by dynamics.** A gapless *non-conformal* theory (`θ≠0`, a non-conserved scalar with power-law decay and no diffusive pole) satisfies "absence of a diffusive mode" yet is gapless, defeating the implication. The "gap detector" is really a **scale detector**: `β(g)≠0` forces the scalar *non-conserved* (CONTINUUM-VARIANCE-DEFECT §4: "an anomaly forces existence, not heaviness"; GAP-IN11 §4), never forces `Γ₀≥λ`.

## Five-point scorecard
1. **Nonperturbative — PASS.** Mori–Zwanzig is an exact IR resummation; all `β`, no loop/`1/N`.
2. **Sign-controlling LOWER bound — FAIL.** Mode *classification*, not an inequality; transport positivity is stability/unitarity, not a bound on `m_eff`; "no diffusive pole" is necessary-not-sufficient for `Γ₀≥λ`.
3. **Non-circular — FAIL.** `Γ₀=G″=m_g`; relaxation-vs-diffusive *is* the gapped/gapless dichotomy.
4. **Non-SUSY + fails free field — PASS (via trace anomaly).** Conformal `T^μ_μ=0` ⇒ conserved ⇒ gapless diffusive ⇒ `G″=0`, cleanly; but the discriminant is conformality (a scale), and a gapless non-conformal CFT fools it.
5. **Uniform `≥λ` — FAIL.** `Γ₀=m_g` is H1; no independent `λ`; conserved-sector diffusion is mass-blind.

## Exact failing estimate
```
Mori–Zwanzig:  ∂_t A = iΩA − ∫K(t−t')A ,  A = Tr F² (0⁺⁺ scalar)
conserved T^{μν}: ω=−iDk²  ⇒ gapless diffusive   (present free-massless / free-massive / gapped YM)
non-conserved scalar: ω=−iΓ₀ , Γ₀=(1/χ)∫₀^∞ dt⟨Ė Ė⟩_disc > 0  ⇒ relaxation
dictionary:  Γ₀ = −b'(0) = G″(k→0) = S−D = 1/Var(a)        [D2; identity]
  "relaxation gap Γ₀≥λ"  ⟺  "G″(k→0)≥λ"  ⟺  Var(a)≤1/λ     ⟺  H1
free massless: T^μ_μ=0 ⇒ conserved ⇒ Γ₀=0 ⇒ G″=0             (acid test PASS)
YM: θ^μ_μ=(β/2g)Tr F²≠0 ⇒ non-conserved ⇒ Γ₀=m_g            (the mass; circular)
missing:  Γ₀≥λ from EFT data alone.  "no diffusive pole" ⇏ Γ₀≥λ;  Γ₀ = pole = gap.
```

**Sharpest line.** The long-wavelength EFT of the `0⁺⁺` sector has *already been forged*: it is the block-spin diffusion with drift `b(a)=G″a` and curvature `b′(0)=−G″` (D2). Hydrodynamics renames its `k→0` relaxation rate "the relaxation gap `Γ₀`" and observes that a gapless theory has none — which is exactly the statement `G″(k→0)>0` = H1. What genuinely *moves* is mass-blind: the conservation-mandated sound/shear modes are present in free and gapped theories alike, and the one discriminating structure — the trace anomaly `θ^μ_μ≠0` — tells you the scalar is non-conserved because a *scale* exists, not that its relaxation rate clears `λ`. Free massless is correctly failed (conformal ⇒ conserved ⇒ diffusive ⇒ `G″=0`), but it is failed by *conformality*, and a gapless non-conformal theory walks through the same door. **#1,#4 pass; the two fails that are #2,#3,#5 are the entire Clay content: a lower bound `Γ₀=G″(k→0)≥λ` does not exist in the hydrodynamic alphabet — it is H1 wearing a fluid's name. NOT 5/5.**

*Refs.* (external) Landau–Lifshitz, *Fluid Mechanics* (1987); Kadanoff–Martin, Ann. Phys. **24** (1963) 419; Forster, *Hydrodynamic Fluctuations, Broken Symmetry, and Correlation Functions* (1975); Mori, Prog. Theor. Phys. **34** (1965) 399; Zwanzig, *Nonequilibrium Statistical Mechanics* (2001); Hohenberg–Halperin, Rev. Mod. Phys. **49** (1977) 435; Son, PRL **98** (2007) 020604 (conformal sound `c_s²=1/3`); Collins–Duncan–Joglekar, PRD **16** (1977) 438 (trace anomaly). (internal) D2-COUPLING-DOEBLIN (`b′(0)=−G″`, `G″=S−D=1/Var`, von Renesse–Sturm); CONTINUUM-VARIANCE-DEFECT §4 ("anomaly forces existence, not heaviness"); GAP-IN11 §4; ABC2-SYNTHESIS (open cell = typical field = H1; `G″(k→0)` target); F2-KREIN-RUTMAN-QUANT / WIENER-HOPF (shape-correct, content-empty; rate `=H1`); AQFT-CAUSAL (mass-blindness; free massless as structural counterexample); C8-SPECIES-NOGO (empty-referent / false-implication pattern); GAP-G2-magnetic (dual confining sector independent of the neutral gap, `W1_0⇏W1_1`). No fabricated bound.
