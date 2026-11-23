# HAMILTONIAN-CHAOS — the transfer operator's quantized chaotic map yields scale-free spectral *statistics* (Wigner–Dyson bulk, soft Airy edge), not the *location* of `E₁`; Gutzwiller is the same `1/ℏ` resurgent sum over the Savvidy saddles — the "classical gap" is `e^{−S_sav}` with a free Maslov sign

> **判词 (verdict).** **No.** Quantum chaos controls *spacing/fluctuation statistics* — scale-free, unfolded, bulk — while the `0⁺⁺` gap `E₁−E₀`, and the narrowed target `G″(k→0)=m_eff²(k→0)≥λ` it encodes, is a *location/scale*. Chaos is **epistemically hostile** here: its hallmark is the *loss of individual-level semiclassical control*, and the edge is exactly where Gutzwiller is invalid (`E→∞`/`ℏ→0`) and Wigner–Dyson universality has not set in. A maximally chaotic system (SYK, a black hole) has an *arbitrarily small* gap: chaos fixes the scrambling *rate*, never the level *spacing*. The one individual-level handle (EBK torus quantization) needs **integrability**; strong-coupling YM is non-integrable. The "classical gap = minimal periodic-orbit action" is `S_sav=8π²/g²`, so the quantum gap is `σ e^{−S_sav}` with `σ` = Maslov/Stokes **free on ℝ⁴** — Gutzwiller *is* resurgence in `1/ℏ` coordinates. Score **2/5**.

## §1 What chaos controls — spacing, not location
`T=e^{−aH}` (`H` = Kogut–Susskind `(g²/2a)ΣE²+(1/g²a)Σ□`) is a positive contraction; `Δ=−lim(1/L)log(λ₁/λ₀)`. Its classical limit is the lattice-gauge flow on `(A,E)` — chaotic (Biro–Matinyan–Müller: positive Lyapunov, KS entropy, Wigner–Dyson bulk). But WD `P(s)~s^βe^{−cs²}` is in *unfolded* units `s=Δ/Δ̄`: it suppresses small spacings — a *normalized, statistical* statement — never `E₁−E₀>λ` absolutely. The mean density is the Weyl volume law, **chaos-independent**; chaos reshapes only its fluctuations, only in the semiclassical window — and the gap lives at the window's bottom, not in it.

## §2 The semiclassical handle dies exactly where the gap is
An *individual* level needs **EBK torus quantization** — actions, `E_n=E(n+α)` — which exists only for **integrable/near-integrable** systems. For chaotic ones EBK fails and **Gutzwiller gives only the *density***, an averaged asymptotic object: YM's non-integrability *removes* the tool that would compute `E₁` from classical data. The harbor (integrable ⇒ EBK ⇒ gap) and the target (strong coupling ⇒ maximal chaos ⇒ no EBK) are **mutually exclusive**; the KAM corner is the weak-field/free-like regime, no gap.

## §3 The "classical gap" is the Savvidy saddle action with a free sign
Gutzwiller `ρ(E)=ρ̄+Σ_pΣ_r(T_p/πħ)|det(M_p^r−1)|^{−1/2}cos(rS_p/ħ−πμ_p r/2)`. The minimal-action orbit in the Savvidy sector is the negative-mode bounce, `S_min~S_sav`, so `Δ~σ e^{−S_sav}`: Savvidy-scale *magnitude*, sign `σ` = Maslov/Stokes **free on ℝ⁴** — verbatim RESURGENCE §2–3. The `1/ℏ` sum is a resurgent trans-series over saddles; the gap identification is its trans-series reading, Stokes sign = H1.

## §4 Chaos is anti-correlated with the gap in the Savvidy sector
The generator of YM chaos **is** the Savvidy negative mode (instability seeds ergodicity) — the *same* `−D²+2ad F` bottom `−g|B|`, count `∝(gB)²V`, both growing with the field (PROOF-SAVVIDY-EXACT-SIGN §1): more chaos ⇔ more negative weight ⇔ **harder** gap. The Ruelle–Pollicott gap is the classical *mixing rate*, matching the *resonance* spectrum only as `ℏ_eff→0` (Faure–Sjöstrand), not the physical gap at fixed `a`; the flow is not uniformly hyperbolic (Savvidy saddle).

## §5 Free-field acid test — passes vacuously
Free massless = free oscillators = **integrable**, Berry–Tabor Poisson, `ω_min=|k|_min→0` gapless. The chaos hypothesis is *false* on it → the tool is silent, no spurious gap (criterion 4 met) — but trivially (integrable), as in RESURGENCE §3 / E4. QE/ETH add nothing: both are *bulk*, density-1, high-`E`, `ℏ→0` statements about *eigenstate* structure, silent on the edge; the edge is a **soft Airy** `ρ~√(E−E₀)` (zero density, no atom, DYSON-GAS §3).

## Five-point scorecard
1. **Nonperturbative — PASS (unbuilt).** Periodic-orbit/instanton actions `O(1/g²)`; genuinely nonperturbative.
2. **Sign LOWER bound `G″(k→0)≥λ` — FAIL.** Statistics/spacing, scale-free & normalized; edge soft Airy, no atom; maximal chaos ⟹ any gap.
3. **Non-circular — FAIL.** Classical gap = Savvidy action, sign = Stokes = H1; RP correspondence only `ℏ_eff→0`; edge orbits need edge data.
4. **Non-SUSY + fails free field — PASS (vacuous).** Non-SUSY ✓; free field integrable/gapless, tool silent ✓ — but never bites at the edge.
5. **Uniform `λ` — FAIL.** WD scale-free; soft edge; no `λ` over `β∈[1/48,2.5]`.

**Exact failing estimate / missing object.**
```
Transfer:  T=e^{−aH}, Δ=−lim(1/L)log(λ₁/λ₀), H = g²/2a ΣE² + 1/g²a Σ□.
Classical: lattice-gauge flow chaotic (λ_L>0, h_KS>0) — controls h_KS / mixing, NOT spec location.
Gutzwiller:ρ(E)=ρ̄+Σ_pΣ_r(T_p/πħ)|det(M_p^r−1)|^{−1/2}cos(rS_p/ħ−πμ_p r/2)
           ρ̄ = Weyl(phase-space vol) = chaos-INDEPENDENT; valid E→∞ / ħ→0 (NOT at E₀).
Edge:      soft Airy ρ~√(E−E₀) ⇒ zero density ⇒ no atom f₁²=0 (DYSON-GAS §3).
WD/bulk:   P(s)~s^βe^{−cs²}, s=Δ/Δ̄ normalized ⇒ suppresses small spacing, NOT absolute Δ>λ.
EBK handle:exists only integrable/near-integrable; YM strong coupling non-integrable ⇒ absent.
Classical gap = S_min ~ S_sav=8π²/g² ⇒ Δ~σe^{−S_sav}, σ=Maslov/Stokes FREE on ℝ⁴.
Max-chaos: SYK/black hole maximally chaotic yet gap arbitrary ⇒ chaos ⇏ gap, no lower bound.
RP gap:    classical mixing rate; =quantum resonance only as ħ_eff→0; not uniform-hyperbolic.
Free field:integrable, Poisson, ω_min→0, gapless ✓ — tool silent (vacuous pass).
```

**Sharpest line.** Quantum chaos is *epistemically hostile* to a gap bound: its one universal theorem (Wigner–Dyson) is scale-free; its one individual-level handle (EBK) needs the integrability strong-coupling YM lacks; its one nonperturbative fragment (Gutzwiller) sums over the very Savvidy saddles whose action fixes a *magnitude* and a *free sign*. Maximal chaos is compatible with *any* gap (SYK), so chaos cannot lower-bound it. The harbor (integrable ⇒ EBK gap) and the target (chaotic strong coupling ⇒ gap) are mutually exclusive: "classical gap vs quantum gap" is not a bound but the trans-series identification `Δ~σe^{−S_sav}` of RESURGENCE, whose Stokes sign is H1. The missing object is an EBK/Gutzwiller handle valid *inside* the chaotic regime — which does not exist — or the stable nonperturbative vacuum with a known Hessian (the atom), which *is* H1.

*Refs.* Gutzwiller, *J. Math. Phys.* **12** (1971) 343; *Chaos in Classical and Quantum Mechanics* (Springer, 1990); Bohigas–Giannoni–Schmit, *PRL* **52** (1984) 1; Berry–Tabor, *Proc. R. Soc. Lond. A* **356** (1977) 375; Shnirelman, *Usp. Mat. Nauk* **29** (1974) 181; Zelditch, *Duke Math. J.* **55** (1987) 919; Colin de Verdière, *CMP* **102** (1985) 497; Deutsch, *PRA* **43** (1991) 2046; Srednicki, *PRE* **50** (1994) 888; Ruelle, *PRL* **56** (1986) 405; Pollicott, *Invent. Math.* **81** (1985) 413; Faure–Sjöstrand, *CMP* **308** (2011) 325; Biro–Matinyan–Müller, *Chaos and Gauge Field Theory* (World Sci., 1994); Savvidy, *PLB* **71** (1977) 133; Nielsen–Olesen, *NPB* **144** (1978) 376; Tracy–Widom, *CMP* **159** (1994) 151; Bleher–Its, *Ann. Math.* **150** (1999) 185. Internal: RESURGENCE; DYSON-GAS; PROOF-SAVVIDY-EXACT-SIGN; C8-SPECIES-NOGO; E4-RMT-EDGE; PROOF-P3-UNCONDITIONAL; PROOF-COARSE-POINCARE (H1); GAP-G2-edge; PROOF-F1SQ; INEQ-R1; W6-BE-at-2.5; CONTINUUM-ASYMPTOTIC-BEYOND; PROOF-SELFDUAL-BACKGROUND.
