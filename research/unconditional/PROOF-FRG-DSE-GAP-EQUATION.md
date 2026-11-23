# PROOF-FRG-DSE-GAP-EQUATION — the gluon gap equation: truncated, gauge-fixed, non-transferable

> **Verdict.** The gap equation is precise but never exact at finite `N`: it closes only under an Ansatz for the three-gluon vertex, so `m²=Π(m²)` can have a rigorous positive solution *only for a specified truncation*, not for Yang–Mills — and never "at β=2.5," where the only defined object is the finite lattice. Existence is not the wall: a model fixed point is standard. The wall is transfer. A rigorous `m_IR>0` is a property of the sign-indefinite, non-reflection-positive Landau-gauge propagator, and it does not enter `H_phys` (the gluon is not in BRST cohomology). **No rigorous physical transfer exists; the obstruction is exactly that of PROOF-LANDAU-GAUGE-IR.**

## 1. The equation, written precisely

Landau gauge, transverse part. The gluon DSE reads `D_T⁻¹(k)=k²+Π_T(k)` with (Roberts–Williams 1994; schematically, full vertices)
```
Π_T(k) = ½ g²C_A ∫d⁴q/(2π)⁴ Γ^{(0)}_{μαβ}D_{αα'}(q)D_{ββ'}(k−q)Γ_{να'β'}
        −  g²C_A ∫d⁴q/(2π)⁴ (…) D_G(q)D_G(k−q) (…)
```
(gluon loop with the *full* `Γ`, ghost loop with the full ghost-gluon vertex; `D_G` the ghost). In the decoupling solution `D_T=1/(k²Z+m²(k²))`, `D_T(0)=1/m²<∞`; imposing the identity at `k=0` gives the **gap equation**
> `m² = Π_T(0)[m²]`,

a nonlinear fixed-point equation (Cornwall 1982; PT–BFM form `m²∝∫y²C(y)Δ²(y)`, Aguilar–Binosi–Papavassiliou 2008; Aguilar–Ferreira–Papavassiliou 2022). FRG is the same object: Wetterich's `k∂_kΓ_k=½STr[(Γ^{(2)}_k+R_k)⁻¹k∂_kR_k]`, and `m²=Γ^{(2)}_A(0)` is the accumulated `k→0` self-energy (Wetterich 1993; Cyrol–Fister–Mitter–Pawlowski–Strodthoff 2016, Pawlowski's vertex-expansion FRG).

## 2. Can a rigorous m>0 be established at β=2.5?

Not for Yang–Mills. `Π_T` is an *infinite* hierarchy: the 2-point sources the full `Γ` (and the 4-gluon sunset), which no Ward identity constrains; no finite-`N` closure exists (the Makeenko–Migdal wall, `cpb-attack/10-dse`). "The gap equation" is thus the equation *after* truncation. For a fixed truncation, existence is real but subtle: a Banach contraction around `m=0` would prove the *trivial* fixed point is unique — the wrong answer; the IR-singular kernel of the Schwinger mechanism (vertex pole) instead makes `Π_T(0)` divergent at `m=0`, so the trivial solution is excluded and a nonzero `m` is forced, a genuine model theorem (Krasnoselskii/Krein–Rutman fixed point of a positive compact operator). But the pole that forces `m>0` *is* the order parameter; assuming it to derive it is circular. And "at β=2.5" is a category error: β is a lattice coupling, the continuum equation carries no β; the only rigorous finite-β object is the unclosed lattice hierarchy, the only rigorous continuum object needs `a→0` uniform control — the crossover certificate. No truncation-independent a-priori bound yields `m>0`.

## 3. The three obstacles

**(a) Gauge dependence.** `D_T` is Landau-gauge, hence gauge-fixing- and Gribov-copy-dependent (Singer: no global section). Its Källén–Lehmann measure is sign-indefinite — Oehme–Zimmermann `∫ρ_A=0` — so a "pole" in `D_T` is not a positive spectral weight.

**(b) Truncation dependence.** The output `m` is the three-gluon vertex's input. The STI fixes only the longitudinal vertex; the transverse vertex carrying the massless Schwinger pole is unprotected (no anomalous slot, `W3-marginal-BPHZ`). Different vertex models give decoupling (`m>0`) or scaling (`m=0`); "controlled" is empirical lattice agreement (Cyrol et al. 2016; Huber 2020; Huber–von Smekal 2013), not proof.

**(c) Transfer — it fails.** Reaching `m_{0⁺⁺}` needs a BRST-closed, positive-norm interpolator with spectral edge at that scale. The gluon is not in `H_phys=ker Q/im Q` (Kugo–Ojima; no colored asymptotic states), and every gauge-fixed composite inherits the sign-indefinite `ρ_A`. No theorem `m_IR>0 ⇒ m_{0⁺⁺}>0` exists; the "glueball ≈ few×m" relation is a lattice coincidence, and Kugo–Ojima's `u(0)=−1` criterion is *violated* by precisely the decoupling branch carrying `m` (`u(0)≈−0.7…−0.8`, Fischer 2008) — the transfer tool is weakest where the mass lives. The failure is thus identical to the Landau-gauge obstruction: gauge-fixed dynamical mass, minus positivity, minus cohomology membership.

**Sharpest line.** The fixed point `m²=Π_T(m²)>0` is real and admits a rigorous *model* existence theorem, but it lives on the gauge-fixed, sign-indefinite, non-cohomological side of the wall; nothing in the equation or in BRST transports it to `m_{0⁺⁺}`, so a rigorous gap-equation solution is evidence for dynamical mass generation — the mechanism — not the gauge-invariant mass gap.

*Refs.* Cornwall PRD 26 (1982) 1453; Roberts–Williams Prog. Part. Nucl. Phys. 33 (1994) 477 (hep-ph/9403224); Aguilar–Binosi–Papavassiliou PRD 78 (2008) 025010 (arXiv:0802.1870); Aguilar–Ferreira–Papavassiliou PRD 105 (2022) 014030 (arXiv:2109.12712); Wetterich PLB 301 (1993) 90; Berges–Tetradis–Wetterich Phys. Rept. 363 (2002) 223 (hep-ph/0005122); Cyrol–Fister–Mitter–Pawlowski–Strodthoff PRD 94 (2016) 054005 (arXiv:1605.01856); Huber Phys. Rept. 879 (2020) 1 (arXiv:1808.05227); Huber–von Smekal JHEP 04 (2013) 149; Kugo–Ojima PTPS 66 (1979) 1; Oehme–Zimmermann PRD 21 (1980) 1661; Singer CMP 60 (1978) 7; Fischer arXiv:0812.2745 (decoupling `u(0)`); Fischer J. Phys. G 32 (2006) R253 (hep-ph/0605173) (review); Maas Phys. Rept. 524 (2013) 203; ATOM-FRG-mass; W6-induced-mass; GAP-G2-mass; PROOF-LANDAU-GAUGE-IR; W3-marginal-BPHZ; cpb-attack/10-dse.
