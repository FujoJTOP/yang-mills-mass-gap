# BROWNIAN-LOOP-SOUP — no positive loop-soup representation of 4D pure YM exists: the exact stochastic-holonomy measure is *2D* (Brownian motion on `G`) or small-λ, coupling a soup to Wilson holonomies makes it **signed** (Weingarten / C8), "connectivity ⇒ gap" is the gap restated, and even a forged mass is `√σ` on the **1-form** sector, not the neutral `G″(k→0)`

> **2026-11-23 · NEW SPECIES (Brownian loop soup / stochastic holonomy).** 判词 (verdict). The Brownian loop soup (BLS) is the Poisson gas of Brownian loops with intensity `c·μ_loop` (Lawler–Werner) — the Symanzik/random-loop representation of a **scalar** field, manifestly *positive*. 4D pure SU(N) YM has **no** positive loop-soup representation: the exact gauge holonomy-loop measure lives only in **2D** (Sengupta, Lévy) and small-λ (never β=2.5), and coupling a soup to Wilson holonomies `w(γ)=(1/N)Re tr W_γ[U]` is **signed** (non-abelian Weingarten) — exactly the positive regrouping PROOF-CONFINEMENT-RANDOM-SURFACE says is missing. "Long loops ⇔ gapless" makes connectivity **the gap restated (#3)**, and the mass a soup could yield is `m₁=√σ`, the flux-tube / **1-form** scale, *not* the `0⁺⁺` curvature `G″(k→0)`. **Acid test PASSES** (free massless = long loops = gapless ⇒ `λ=0`). **#1 nominal, #4 PASS; #2,#3,#5 FAIL = 2/5. NOT 5/5.**

## §1 The forge and the dictionary
Target: `G″(k→0)≥λ`, with `G″=S−D=1/Var(a)` the `0⁺⁺` coarse-mode curvature and `G″(k→0)≥λ ⟺ Var(a)≤1/λ` = **H1** (OTHER-H1-POINCARE, D2). The BLS loop of duration `t` carries a killing `e^{−m²t}` *only when the scalar has mass m*; its length is exponentially cut off iff `m>0`. Two YM couplings: (i) **stochastic holonomy** `Pexp∮_γ A` of `U~μ_YM`; (ii) **gauge-dressed soup** weighting loops by `w(γ)`. Question: does the loop-density / holonomy-distribution decay give `m_eff≥√λ`?

## §2 Where the exact object lives — not 4D β=2.5
The Symanzik loop gas is **scalar** (Symanzik 1969; Le Jan); the gauge measure expands over **surfaces with signed Weingarten weights** (CPS; PROOF-CONFINEMENT-RANDOM-SURFACE) — and that signed surface sum is what the *positive* soup would have to be. The one exact gauge holonomy-loop measure — stochastic holonomy = Brownian motion on `G`, heat-kernel `⟨W(C)⟩=p_{area}(e)` (Sengupta Mem. AMS **126**/600; Lévy Mem. AMS **166**/790; Driver) — is **2D**: it proves the 2D area law, but 2D YM has no 4D gap to transport. In 4D at β=2.5 the holonomy decay **is** the area law = the open large-λ problem (CNS: λ<λ₀=1/24 only; β=2.5 sits ~75× above). A coupling (i)/(ii) makes `w(γ)` a *signed* non-abelian character; percolation/connectivity (FKG, Harris, DPP) needs a **positive** measure — absent on torsion SU(N) (C8). The soup "density" is not a probability.

## §3 Connectivity is the gap, and the wrong sector
Gapped ⟺ loops have an exponential length cutoff ⟺ no large-scale percolation ⟺ the two-point function decays. So "long loops ⇔ gapless" **is** the mass gap; the decay rate to extract *is* `m`, readable only from the measure = the field = the gap. **Circular (#3).** Worse (the sharpest miss): a soup probes **lines**, so its IR scale is the long Wilson-loop → string tension `σ` → flux-tube `m₁=√σ`. The target `G″(k→0)` is the **`0⁺⁺` plaquette/surface (two-form)** coarse curvature, neutral sector; the area law does **not** imply the `0⁺⁺` gap (`m_{0⁺⁺}=c√σ` measured, not theorem — PROOF-CONFINEMENT-RANDOM-SURFACE, PROOF-STRING-TENSION-STABILIZER). The soup's "mass" is on the **wrong sector** (1-form/center), blind to `tr F²` — the same cross-channel immunity as dual-Meissner / the Polyakov holonomy (k0-HOLONOMY: 1-form, free-gluon `V″∝1/L²`).

## §4 Free-field acid test (passes, on the free side)
Free massless: `μ_loop` has `p_t(0,0)~t^{−d/2}`, no `e^{−m²t}` ⇒ **long loops, no exponential length decay** (d=4 marginal `∼log`), no large-scale shutoff ⇒ gapless ⇒ `λ=0`. The tool correctly returns **no gap**: **acid test PASSED** (non-SUSY, no fabricated mass — #4). Massive free: `e^{−m²t}` cutoff, "mass" = the *input* `m`. That is the whole content: the soup reads the scalar's **input** mass, never an emergent YM gap — and it passes #4 *because it is the free field's own representation*; there is no positive YM soup to be nonzero on.

## Five-point scorecard
1. **Nonperturbative — PASS (nominal).** Soup and 2D stochastic holonomy are exact/nonperturbative; the 4D β=2.5 representation does not exist (2D / small-λ only).
2. **Sign-controlling `G″≥λ` — FAIL.** Soup positive, YM signed (Weingarten/C8); no positive domination ⇒ no lower bound.
3. **Non-circular — FAIL.** Connectivity ⇔ long-loop ⇔ gapless ⇔ mass is definitional; the mass is `√σ` (1-form), not `G″(k→0)` (0⁺⁺).
4. **Non-SUSY + fails free field — PASS (formal).** Free massless long-looped gapless ⇒ `λ=0`; but it is the free object, no YM arm.
5. **Uniform `≥λ` — FAIL.** 2D exact does not transport; 4D soup absent at β=2.5; percolation threshold model/dimension-specific; no uniform λ.

## Exact failing estimate
```
Target:  G″(k→0) ≥ λ   ⟺   H1 (Var(a) ≤ 1/λ).
BLS: μ_loop = ∫₀^∞ dt/t ∫dx P^t_{x,x}(dγ);  mass m ⇒ length ~ e^{−m²t}  (input mass only).
  • scalar/Symanzik object; 4D gauge YM = SIGNED Weingarten surfaces (CPS); no positive soup.
  • exact gauge holonomy-loop measure: ONLY 2D (Sengupta/Lévy, Brownian motion on G, heat kernel).
  • 4D holonomy decay ⇔ area law = open large-λ (CNS λ<1/24; β=2.5 ~75× above).
  • coupling w(γ) ⇒ SIGNED (C8: no positive Gibbs measure on torsion SU(N)); density not a probability.
connectivity: long-loop ⇔ gapless ⇔ gap                       [definitional, #3]
sector:      soup probes 1-form lines ⇒ m₁=√σ ≠ G″(k→0) (0⁺⁺, neutral).   [wrong sector]
free massless: no e^{−m²t} ⇒ long loops ⇒ λ=0.              [ACID PASS, free side]
missing:  a POSITIVE 4D YM loop measure at β=2.5 (does not exist) + 0⁺⁺ transport from a 1-form mass.
```

**Sharpest line.** The loop soup is exactly the **positive measure** PROOF-CONFINEMENT-RANDOM-SURFACE says is missing; it exists for the free field *because* the free field is positive, and dies on YM at the non-abelian Weingarten sign (C8). "Long loops ⇔ gapless" turns connectivity into the gap restated — and even a forged decay delivers `√σ` on the **1-form** sector, not the neutral `0⁺⁺` `G″(k→0)`. The genuine exact object (stochastic holonomy, Sengupta/Lévy) is 2D; its 4D analogue is the open large-λ area law, not a lever. **#1 nominal, #4 formal; #2, #3, #5 — the Clay content — fail. No fabrication.**

*Refs.* (external) Lawler–Werner, *Probab. Theory Related Fields* **128** (2004) 565–588 (BLS); Lawler–Schramm–Werner, *JAMS* **16** (2003) 917; Symanzik, in *Local Quantum Theory* (1969); Le Jan, *Ann. Probab.* **38** (2010) 1280; Sengupta, *Gauge Theory on Compact Surfaces*, Mem. AMS **126** (1997) no. 600; Lévy, *Yang–Mills Measure on Compact Surfaces*, Mem. AMS **166** (2003) no. 790; Driver, *J. Funct. Anal.* **110** (1992) 272; Migdal, *Sov. Phys. JETP* **42** (1976) 413; Makeenko–Migdal, *Phys. Lett. B* **88** (1979) 135; Cao–Park–Sheffield arXiv:2307.06790; Cao–Nissim–Sheffield arXiv:2505.16585, arXiv:2509.04688. (internal) PROOF-CONFINEMENT-RANDOM-SURFACE; E3-LOOPSPACE; C8-SPECIES-NOGO; k0-HOLONOMY-EFFPOT; PROOF-ROUGH-PATH-SIGNATURE; OTHER-H1-POINCARE; D2-COUPLING-DOEBLIN; INVENT-VARIANCE-INEQ; PROOF-STRING-TENSION-STABILIZER.
