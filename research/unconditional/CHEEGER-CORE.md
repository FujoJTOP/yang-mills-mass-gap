# CHEEGER-CORE — 9th-species CORE attack: Cheeger's isoperimetric spectral-gap inequality on the YM measure. **STEP-0 finding: Cheeger/isoperimetry is already banked** (`B3-ISOPERIMETRIC` 2/5; classified a re-tread in `OUT-transplant` §11; Dirichlet-form-conductance sibling `POTENTIAL-THEORY-CAPACITY` 2/5; RCD sibling `RCD-SYNTHETIC-RICCI` 2/5). Cheeger is **NOT a new species**. The one genuinely-untested question — *does the "negative curvature can still give h>0 (hyperbolic/expander)" escape survive the **saddle-type** Savvidy negativity?* — is decided **NO**: expanders buy `h>0` from volume **expansion**; the Savvidy `−g|B|` is an **unstable (saddle) neck** that *spreads* the measure, so isoperimetry **sees the bottleneck more directly than Bakry–Émery did**. **Conforms: 2/5.**

> **2026-11-23 · 判词 (verdict). NO NEW SPECIES; wall certified.** For the Langevin/OU generator of `μ=Z^{−1}e^{−S_W}dU` on the gauge quotient `M=G^E/G_site`, `S_W=(1/2g²)∫trF²`, Cheeger gives the **curvature-free** gap bound `λ₁ ≥ h²/4`, `h=inf_A μ⁺(∂A)/min(μ(A),μ(Aᶜ))`. This is the campaign's best structural hope — it is **not** `Γ₂/ Ric`, so it does not literally die at the Savvidy point `min(Ric+Hess S)≈−2.3` where Bakry–Émery died. **It dies one step later, at the same object.** Lower-bounding `h` is, **in content**, bounding the transverse spread of the soft `k→0` mode: `h ≥ h₀ uniform ⟺ Var_μ(a_{k→0}) ≤ 1/λ uniform = H1`. The precise geometry: expanders (hyperbolic quotients) have `h>0` because their negative curvature is **volume-expanding** (balls grow exponentially, no thin neck); the Savvidy `−g|B|` (min `≈−2.3` at `β=2.5`, 7.1% transverse-negative, `W4`) is a **saddle** — `S` concave in the transverse-magnetic direction ⇒ the measure is **widest** there ⇒ the cheapest halving set sits on that neck ⇒ `h_binding ∼ λ^{1/2}` and Cheeger becomes `λ ≥ (λ^{1/2})²/4 = λ/4`, a tautology. **The exponential-neck hope is nonetheless refuted as "hard":** the Arrhenius barrier at the target shell `|a|≤λ^{−1/2}` is `H∼0.3/λ=O(1)` (POTENTIAL-THEORY-CAPACITY §3), so the profile does **not** collapse exponentially — the failure is a **soft bootstrap on `λ` = H1**, not a wall. Passes #1, #4 (clean, non-vacuous); fails #2, #3, #5. **2/5 — identical to B3/CAPACITY/RCD; a re-tread, not a 9th species.**

---

## §0. STEP-0 honesty discharge (what is already banked)

The corpus already contains the full Cheeger/isoperimetric analysis and its two siblings; this report does not re-derive them.

- **`B3-ISOPERIMETRIC.md`** (2026-11-23): Cheeger on the **coarse marginal** `ν=P_*μ`, `ρ(ν)≥h(ν)²/4`. Verdict "no constructive ρ": the **two-coordinate obstruction** — `μ=ν(da)·μ(·|a)`, marginal and conditionals are independent coordinates, so no fiber property bounds `ρ(ν)`; `ν`-Cheeger competitors are **fiber-saturated** (`P̄^{−1}(B)`), whose fiber-tangent boundary is `≡0`, so the fiber LSI is *silent* on exactly the sets that define `h(ν)`. Residual `h(ν)≥2√ρ ⟺ transverse bottleneck bounded ⟺ H1`. Score **2/5**.
- **`POTENTIAL-THEORY-CAPACITY.md`** (2026-11-23): the capacity/conductance form `λ₁≥Φ_eq²/4` (Cheeger in capacity dress, `λ₁≥h²/4` with `h=Φ_eq`). Cleanest **non-vacuous** free-field discriminator (flat zero mode ⇒ `Cap=0`). Score **2/5**; uniform floor `Φ_eq≥c` = H1.
- **`RCD-SYNTHETIC-RICCI.md`** (2026-11-23): Lott–Villani–Sturm `RCD(K,∞)`; on flat `A` synthetic = pointwise `Hess`, `K_*=−g|B|+N/12<0`; compactness-Poincaré `ρ∼1/D²→0`. Score **2/5**.
- **`OUT-transplant.md` §11** (prior art): already **classifies Cheeger as a re-tread** — "Cheeger/expander bounds are conductance aggregates of the same kind. Re-tread of P4/T1/T4."

**Therefore:** "Cheeger = 9th species" is the **isoperimetric/conductance/RCD cluster** re-entered. Per project rule, do **not** re-derive. What follows isolates the one claim none of the four makes: the **saddle-vs-expander** refutation of the specific hope stated for this species (isoperimetry escapes where BE died), and the **full-`μ` vs marginal-`ν`** bookkeeping. Verdict is graded, not inflated.

---

## §1. Cheeger / Buser for the YM generator, with `h` written out

**Generator and measure.** `L = Δ_M − ∇S_W·∇_M` on `M=G^E/G_site` with bi-invariant metric, invariant `μ=Z^{−1}e^{−S_W}dU`, `S_W=(1/2g²)∫trF²` (Langevin/OU; equivalently the transfer operator's Dirichlet form `ℰ(f)=∫|∇f|²dμ`). Spectral gap `λ₁ = inf_f ℰ(f)/Var_μ(f)`.

**Cheeger's inequality** (Cheeger 1970; for general metric-measure spaces, **no curvature hypothesis**):
```
λ₁ ≥ h(μ)²/4,     h(μ) = inf_{A⊂M} μ⁺(∂A) / min(μ(A), μ(Aᶜ)),
```
`μ⁺(∂A)` the Minkowski/content boundary measure (coarea `∫|∇1_A|dμ`). (The task's `h²/2` is the `Δ`-vs-`2Δ` normalization variant; we use the conservative `h²/4` matching B3.) The same constant on the coarse marginal `ν` is `h(ν)`, and on `ν` this reads `ρ(ν)=λ₁(ν)≥h(ν)²/4` (B3).

**Buser's converse** (Buser 1982): for a Riemannian `n`-manifold with `Ric ≥ −(n−1)δ`,
```
λ₁ ≤ C·(h + √h)         (δ=0: λ₁ ≤ C h;  general: λ₁ ≤ 2(n−1)δ + C(h+√h)).
```
**Crucial asymmetry (used in §3):** Cheeger's *lower* bound is **universal** (any space); Buser's *upper* bound needs a **two-sided/lower**-Ricci control. On the YM quotient the relevant curvature is `Ric_Haar+Hess S_W⪰−g|B|+N/12 = K_*<0` **but unbounded below** across the Savvidy cone (indefinite). Hence **Buser does not apply**, so `h` and `λ₁` are **not tied in both directions** — the only usable arrow is Cheeger `λ₁≥h²/4`. This matters: it means one *cannot* read `h>0` off `λ₁>0`; the only route to `h>0` is a **direct perimeter bound**.

**The object is the coarse bottleneck.** `h(μ)` is minimized over **all** sets; the lifted coarse (`ν`-)sets are a sub-family, so `h(μ) ≤ h(ν)` (§2). The binding minimizer is the set separating the soft `k→0` amplitude, whose Cheeger ratio is the transverse spread `∼1/σ`, `σ²=Var_μ(a_{k→0})=1/λ` (D3). By Cheeger/Buser on the **good** directions this is `h(ν) ≍ √ρ`; on the **bad** transverse direction it is the target itself.

---

## §2. Isoperimetry ≠ curvature — can Cheeger escape the Savvidy point? (the central hope, decided NO)

**The hope is half-correct, and we credit it precisely.** Cheeger's inequality carries **no** `Ric`/`Γ₂` hypothesis. General theory: **hyperbolic** `n`-manifolds (`Ric=−(n−1)`) have `h>0`, and **expanders** have `h≫0` with arbitrarily negative curvature. So isoperimetry is genuinely **not** the same gate as Bakry–Émery — BE died *literally* at `Γ₂=Ric+Hess S≈−2.3`; Cheeger does not read that number. **At the level of hypotheses, the escape is real.** This is the one respect in which the 9th-species proposal is sharper than a pure re-tread.

**Why it nonetheless fails for YM — saddle-type, not expander-type, negativity (the new point).** There are **two ways** to have `Ric<0` and `h>0`:
1. **Expander/hyperbolic type (stable-negative):** negative curvature is uniform in *all* directions; ball volume grows exponentially `V(B_r)∼e^{cr}`; every halving cut has perimeter `Ω(1)` per unit measure; **no thin neck**. This is what powers the "negative curvature can still give `h>0`" precedent.
2. **Saddle type (unstable-negative):** curvature is negative in **one** (few) transverse-magnetic direction and positive elsewhere; near the saddle `S≈S_0−|λ_neg|x²/2 ⇒ e^{−S}` is **widest** along that axis. The measure is *spread* in that direction; a halving set placed at the saddle cuts at perimeter `∼` density `× σ`, which **vanishes** as `σ→∞` (`λ→0`). A saddle **is** a thin neck for the measure.

**The Savvidy `−g|B|` is saddle-type, decisively.** `Hess S_W=−D²+2adF_B` indefinite; at `β=2.5` the transverse-magnetic mode is negative (`min≈−2.3` SU2 L=4 / `−2.0` SU3 L=3; `≈−1.1β`; `λ_min=−g|B|`), **7.1% gauge-projected transverse-negative at the SU(2) crossover after quotient** (W4 — real physics, not gauge). These are **literal saddle necks**: the measure fattens along exactly the directions that define the cheap Cheeger competitors. So the hyperbolic/expander precedent **does not transfer**: YM's negativity is not a uniform expansion, it is a set of **7% unstable necks at the `k→0` transverse-magnetic mode**, and isoperimetry — being a *bottleneck functional by definition* — registers them **more directly** than `Γ₂`, which at least had the positive Haar term `N/12` as partial cover. **The escape fails, and fails for a geometric reason sharper than "it's curvature": the specific negative curvature of Yang–Mills is the one kind isoperimetry cannot exploit.**

**Full-`μ` does not rescue it either (bookkeeping).** `h(μ)≤h(ν)`: adding the fiber (gauge-orbit) cut directions can only *lower* the infimum. The gauge direction is the **good** one — compact SU(N), log-concave Haar, `Φ″⪰N/12`, `W_geo≥0` (Prékopa), `β`-independent positive — so it does **not** supply a cheaper cut. The cheapest cut is therefore **transverse** = the `ν`-bottleneck = the same `H1` object B3 already isolated. Going from marginal `ν` to full `μ` **inherits and worsens** the B3 obstruction; it adds no new positive direction.

---

## §3. The two acid tests (both run decisively)

**ACID TEST 1 — does Cheeger give the gap *without* the measure-positivity circularity? Formal NO/YES split (this is the crux).**

- **Non-circular in FORM — YES.** `λ₁≥h²/4` derives the gap from a **set-measure perimeter** `h`, *not* from `Γ₂`/`Ric`/`G″` and *not* from assuming a gap. It is genuinely a different syntactic route than Bakry–Émery.
- **Circular in CONTENT — YES (the kill).** The only available handle on `h` is the spread of the soft mode. For the binding transverse cut with effective width `σ=√Var=λ^{−1/2}`, `h_binding∼1/σ=√λ`. Substituting into Cheeger: `λ₁ ≥ (√λ)²/4 = λ/4` — a **tautology** that yields no positivity. There is **no** independent perimeter lower bound: every direct route to `h≥h₀` reads the same `Var(a_{k→0})` (§2), and Buser cannot convert `λ>0⟹h>0` precisely because it is inapplicable (Ric unbounded below, §1). Formally:
  ```
  h ≥ h₀ (uniform)  ⟺  transverse bottleneck bounded  ⟺  Var_μ(a_{k→0}) ≤ 1/λ (uniform)  =  H1.
  ```
  This is **exactly B3's residual** — derived here from the saddle geometry rather than the two-coordinate independence. **#3: non-circular in form, circular in content ⇒ FAIL at closing.**

**ACID TEST 2 — free field (massless must give `h=0`; massive `h>0`, gap `=m`). CLEAN PASS (non-vacuous — Cheeger's genuine win over BE).**

- **Free massless** (continuum / thermodynamic limit): Gaussian `N(0,C)`, `C=g²(−Δ)^{−1}`; `λ_max(C)=g²/λ_min(−Δ)`, `λ_min(−Δ)∼(2π/L)²→0`. For a Gaussian, `h=√(2/π)/σ_max`, `σ_max=√λ_max(C)⇒ h ∼ √λ_min(−Δ)/g → 0`. **Cheeger returns `h=0`, gap `0` — the gapless theory is correctly failed, non-vacuously** (not "0 trivially from compactness": this is the *continuum* free field, no compactness). This is strictly stronger than B3's fiber-LSI (which was acid-test-*blind*) and matches `POTENTIAL-THEORY-CAPACITY`'s `Cap=0`; the two are the same discriminator (`h=0 ⟺ Cap=0`).
- **Free massive:** `C=g²(−Δ+m²)^{−1}`, `σ_max=g/m`, `h∼(√(2/π))·m/g>0`, gap `=m²/g²`; `h∼m` with `gap=m²` — **correct scaling.** Cheeger tracks the free gap correctly on both sides.

**Verdict on acid tests.** #4 PASSES **cleanly and non-vacuously** (a real advantage over BE/W6-BE, whose "Haar `N/12`" was a compactness artifact leaking on the massive side; Cheeger has no such leak — it returns the *right* `0` for massless and the *right* `m` for massive, in the continuum). But the free-field correctness is **not** a discriminator of the YM **atom**: Cheeger returns whatever gap the measure carries and is **blind to `f₁²`** (a spectral gap is necessary, not sufficient — MS-BLIND §3, SIGNED-FOURIER-BLIND §0). Passing #4 does not touch #2/#3/#5.

---

## §4. The Savvidy-crossover isoperimetric profile — does it collapse?

**Finite `L`:** trivial — `M` compact, `e^{−S_W}>0` smooth ⇒ `h(μ)>0`, `λ₁>0` for every `L`. This is the CAPACITY/B3 "sign is free but vacuous in finite volume" point; all content is `L`-uniformity (#5).

**Uniform-in-`L` profile across the cone — the honest shape.** Along the transverse-magnetic `k→0` axis the local density is `e^{−S}`, `S=S_0−|λ_neg|x²/2+O(x⁴)` (saddle; `|λ_neg|≈g|B|`, min `−2.3` at `β=2.5`), globally coerced by the quartic growth `S∼β‖F‖²V` at large amplitude. Two regimes:
- **At the target thermal shell** `|a|≤λ^{−1/2}`: the Arrhenius crossing cost is `H∼|V̄″|·λ^{−1}≈0.3/λ = O(1)` for `λ=O(1)` (POTENTIAL-THEORY-CAPACITY §3, same estimate). So there is **no exponential thin neck** at the scale where the mode lives; the profile is **controlled up to an `O(1)` barrier** — it does **not** collapse `∼e^{−cL}`.
- **Uniform floor:** the only way to make "no exponential neck" quantitative is to seed `λ≥c`, i.e. `h≥h₀`. That seed **is** H1. So the profile is controlled **conditionally on the bootstrap `λ≥c`** — the same soft bootstrap the capacity route found.

**Answer:** the isoperimetric profile at the Savvidy crossover **does not catastrophically collapse** (no `e^{−H}` neck with `H` growing — the Savvidy dip is `O(1)` at the target shell). It is killed **softly**, by the **bootstrap-on-`λ` seed = H1**, not hard, by an exponential wall. This is the honest "residue": a *milder-looking* (barrier `O(1)`) but *still unclosed* estimate.

---

## §5. Scorecard (0–5), sharpest line, exact stuck step

**Five-point scorecard**
1. **Nonperturbative — PASS.** Cheeger's inequality is an exact measure-geometry theorem, all `β`, no loops; `h` is a global set-measure functional. (Genuine, curvature-free — the structural reason the species looked promising.)
2. **Sign-controlling LOWER bound `m_eff²(k→0)≥λ` — FAIL at the bottleneck.** Cheeger *is* a lower bound (`λ₁≥h²/4`) — the right *shape*, unlike the positivity tools — but bounding `h` from below reduces (saddle, §2) to the transverse spread `=1/λ`; no independent perimeter floor exists. `h_binding∼√λ` gives the tautology `λ≥λ/4`.
3. **Non-circular — FAIL at closing (the precise crux).** Non-circular in **form** (gap from perimeter, not `Γ₂`/`Ric`/gap); **circular in content** (`h≥h₀⟺H1`; Buser inapplicable to break the tie; full-`μ` only adds the good gauge direction, §2). A real syntactic advance over BE, zero substantive gain.
4. **Non-SUSY + fails free massless field — PASS (clean, NON-vacuous).** Pure YM is non-SUSY. Free massless: Gaussian `h∼√λ_min→0`, gap `0` (correctly failed, in the *continuum*, not via compactness); free massive: `h∼m`, gap `m²` (correct). **Genuinely discriminating on the free axis** — sharper than BE's compactness-artifact "pass."
5. **Uniform `(L,a)`, reach `β≈2.5` — FAIL.** Uniform `h₀` is exactly H1; finite-`L` positivity is vacuous; across the cone the profile is controlled only *conditionally on the bootstrap `λ≥c`* (§4), which is H1.

**Score: 2/5** (passes #1, #4; fails #2, #3, #5) — **identical to `B3-ISOPERIMETRIC` (2/5), `POTENTIAL-THEORY-CAPACITY` (2/5), `RCD-SYNTHETIC-RICCI` (2/5).** Cheeger is a **re-tread** of the banked isoperimetric/conductance/RCD cluster, not a 9th species. It is the **best-motivated** member (curvature-free, non-vacuous free-field pass) and still lands where all of them land.

**Sharpest line.** *Cheeger escapes Bakry–Émery's **hypothesis** (no `Ric`, no `Γ₂<0`) but not its **conclusion**: hyperbolicity/expanders buy `h>0` from volume **expansion**, whereas the Savvidy `−g|B|` (min `−2.3`, 7.1% transverse-negative) is an **unstable saddle** that *spreads* the measure — so isoperimetry, a bottleneck functional, sees the soft-mode neck **more directly** than `Γ₂` did. Lower-bounding `h` is, in content, `h_binding∼√λ`, and Cheeger becomes the tautology `λ≥λ/4`. The exponential-neck hope is refuted as *hard* (barrier `O(1)` at the target shell, not `e^{−cL}`); what remains is the *soft* bootstrap `λ≥c` — the seed that **is** H1. Isoperimetry does not cross the wall; it labels it, on the same square the other eight species died on, in the cleanest possible font.*

**Exact stuck step.**
```
target:   λ₁ ≥ λ₀>0  uniform (L,a),  β up to 2.5           [= Poincaré gap ⇒ H1]
Cheeger:  λ₁ ≥ h²/4,   h(μ)=inf_A μ⁺(∂A)/min(μ(A),μ(Aᶜ))    [curvature-free — genuine]
binding:  cheapest cut = transverse k→0 saddle; h_binding ∼ 1/σ = √λ
          ⇒  λ ≥ (√λ)²/4 = λ/4  (tautology — no positivity)
missing:  a perimeter lower bound  h(μ) ≥ h₀  that does NOT read Var(a_{k→0})
          ⟺  h ≥ h₀  ⟺  Var_μ(a_{k→0}) ≤ 1/λ  ⟺  m_eff²(k→0) ≥ λ  =  H1.
          (Buser cannot supply λ>0⟹h>0: needs Ric≥−K, false — Hess unbounded below.)
residue:  §4 profile does not collapse exponentially (H∼0.3/λ=O(1) at target shell);
          it needs only the soft seed λ≥c — which is H1. A soft bootstrap, not a wall.
```

**Refs.** *Internal (banked, cited as ground, not invented):* **B3-ISOPERIMETRIC** (Cheeger on marginal `ν`, `ρ≥h²/4`, two-coordinate obstruction, fiber-saturated `∂≡0`, `h(ν)≥2√ρ⟺H1`, 2/5); **POTENTIAL-THEORY-CAPACITY** (`λ₁≥Φ_eq²/4`="Cheeger in capacity", free `Cap=0`, Arrhenius `H∼0.3/λ=O(1)` at target, bootstrap-on-λ, 2/5); **RCD-SYNTHETIC-RICCI** (`K_*=−g|B|+N/12<0`, compactness `ρ∼1/D²→0`, Buser needs Ric, 2/5); **OUT-transplant §11** ("Cheeger/expander = conductance aggregate, re-tread of P4/T1/T4"); TRANSPORT-BLIND / W6-BE-at-2.5 (`min(Ric+Hess S)≈−2.3/−2.0`, `−1.1β`, `m<0`); B3-FUNCTIONAL-BE (`Γ₂=G″` lossless = H1; `N/12` compactness caveat); PROOF-SAVVIDY-EXACT-SIGN (`−D²+2adF`, `λ_min=−g|B|`, negativity non-BPS); W4-blocked-positivity (7.1% gauge-projected transverse-negative at SU2 crossover); MARGINAL (`PᵀHP` min `−0.3`); OTHER-H1-POINCARE (Haar half, `Φ″⪰N/12`, log-cut anti-correlation); D3-SUSCEPTIBILITY-VARIATIONAL (`G″=1/Var(a)`); PROOF-COARSE-POINCARE (exact coarse form, Lemma SC β<1/48); PROOF-JACOBIAN-SAVES (`G″=Δ_dyn+W_geo+X`, `W_geo≥0` Prékopa); MS-BLIND §3, SIGNED-FOURIER-BLIND §0 (spectral gap ≠ atom, `f₁²`); docs/YANG-MILLS-MASTER.md §5–6 (residual object = variance upper bound, one species at UV/IR). *External (standard):* Cheeger, *Problems in Analysis* (1970) 195; Buser, *Ann. Sci. ENS* **15** (1982) 213; Lawler–Sokal, *Trans. AMS* **309** (1988) 557 (λ≥h²/4 normalization); Buser / Ledoux *The Concentration of Measure Phenomenon* (2001); Gromov–Lévy; Savvidy, *Phys. Lett. B* **71** (1977) 133.
