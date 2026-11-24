# SIGNED-FOURIER-BUILD — the representation-ring / signed-Fourier (Bochner) species, built as a precise inequality: the atom is a non-negative representation coefficient, its LOWER bound is a signed 6j-recoupling dominance statement, and the dominance fails only at the crossover where it is exactly the constituent mass gap `m_g>0`

> 2026-11-23 · 4D pure SU(N) YM · Constructive forge of the **representation-ring / signed-Fourier (Bochner-type)** species for `f₁²=|⟨0|trF²|G⟩|²`, *built as an inequality, not prose*.
> **判词 (verdict).** The construction is real and precise but **NOT 5/5 — 2/5 (1/5 on the strict red-team read; see §4).** The connected `trF²` correlator is a sum over representation-labelled surfaces; its Källén–Lehmann measure is a sum of **non-negative representation coefficients** `w_r≥0` (`w_r=|A_r|²`), i.e. the "signed-Fourier" transform of the connected correlator is a *positive* measure on the representation ring (Bochner/Godement positivity = RP). The atom is the lightest coefficient, `f₁²=w₁`. Promoting `w₁≥0` to `w₁≥r₀>0` is a **dominance statement**: the minimal fundamental tube amplitude must beat the **signed 6j/Weingarten recoupling tail** `𝒯(u)`. The precise inequality is `f₁² ≥ c₁(β)u^p(1−𝒯(u))²₊`. It **(i) reproduces P3 at `β<1/48`** (`𝒯=O(u)` ⟹ `f₁²≥r₀>0`, consistency); **(ii) correctly vanishes on the free field** (the free two-gluon operator has `A₁=0`, the cut edge-weight `→0` as `t^{-3/2}` — acid test PASSED); **(iii) breaks at the crossover `β≈2.5`**, where `u=c_f(2.5)≈2.0>1` (measured) leaves the small-activity radius entirely and `𝒯(u)=O(1)` is uncontrolled — restoring `𝒯(u)<1` is *literally* "the minimal fundamental tube is a non-cancelling isolated bound state" = isolation `I2` = **`m_g>0` = H1**. The 6j symbols are the only signed objects; on `SU(N≥3)` they are sign-indefinite with no termwise sign, so no representation-theoretic positivity survives to `2.5`. Sharpest: the species isolates the atom floor to one named scalar — the signed recoupling tail — whose control is the mass gap, not a representation datum.

---
## §0 The one-line inequality

```
f₁²  =  |A₁|²  ≥  c₁(β)·u^p · ( 1 − 𝒯(u) )²₊ ,        u := c_f(β) ∈ (0,1) at weak activity,
   𝒯(u) := Σ_{γ≠min} |Wg_γ / Wg_min| · u^{|γ|−|min|}    (signed 6j/Weingarten recoupling tail, |6j|≤1).
Hypothesis (H_coh):  𝒯(u) ≤ 1 − ρ  (ρ>0)   ⟹   f₁² ≥ ρ² c₁ u^p =: r₀ > 0.
```
`c₁(β)>0` is the operator-content constant (`f⊗f̄→1` recoupling) and `p` the content exponent (P3's bookkeeping). (H_coh) is the whole game; it is the *precise* statement of "coherence, not disorder."

---
## §1 — Character expansion of `C(t)=⟨O(t)O(0)⟩_c`, `O=trF²`, and the representation coefficient carrying the atom

**1.1 The two expansions of one measure.** On `T⁴_L`, `dU_e` Haar, Wilson action `S=βΣ_p(1−N⁻¹Re tr U_p)`. Two facts fix everything:
- (Itzykson positivity) `e^{(β/N)Re tr U_p}=Σ_r c_r(β)χ_r(U_p)`, with **`c_r(β)>0` for all `β>0`** (Drouffe–Itzykson). `c_0=1+O(β)`, `c_r≍β^{|r|}` at small `β`. *(Verified numerically, SU(2): `c_0=2I₁(β)/β`, all `c_j>0`; `w7/verify_signedfourier.py`.)*
- (RP / Osterwalder–Seiler) For the gauge-invariant `O`, the connected correlator is a positive spectral integral,
  `C(t)=∫₀^∞ e^{−st} ν(ds)`, `ν≥0`, `m₁²=inf supp ν`, **`f₁²=ν({m₁²})`**.

The representation-ring and the spectral measure are the *same object in two bases*. Grouping the connected spin-foam (dual-complex) surfaces `γ` that interpolate the two plaquettes by the representation content they carry, `O=Σ_r ô_r χ_r` (`ô_f=ô_f̄=1/(2N)`, only `f,f̄` at single-plaquette level),

```
C(t) = Σ_r |A_r|²  e^{−m_r t}  +  ∫ e^{−st} ρ_c(s) ds
     ║  (Bochner / signed-Fourier)  ║
ν    = Σ_r w_r δ_{m_r²}  +  ρ_c ds ,   w_r := |A_r|² ≥ 0 ,   A_r = Σ_{surfaces in rep r} (Π character coeffs)·(Π recoupling).
```

**1.2 Which coefficient carries the atom.** The lightest non-vacuum representation is the **fundamental `f`** (minimal nonzero Casimir `C₂(f)=(N²−1)/2N`; verified SU(2) `C₂(½)=¾<C₂(1)=2<…`), so the lightest tube is the single fundamental flux loop and `m₁=m_f`; the atom is `f₁²=w₁=|A_f|²`. This is P3's "one fundamental tube." The amplitude is a **modulus-square of a recoupling network**; `w_f≥0` is Bochner/Godement positivity on the representation ring — it *is* the KL positivity `ν≥0`, nothing more. The "signed" content is only in the intermediate computation: the **6j/Racah/Weingarten symbols** in `A_f`, which alternate.

**1.3 The minimal surface is positive; the tail is signed.** The minimal connected surface is the single `1×t` fundamental tube; closing `f⊗f̄→1` at each end has coefficient `+1` (P3 §3), so the bare tube amplitude is `A_f^{(0)}=a₀ u^{t/a}`, `a₀=c₁u^p>0`, `m_f=(1/a)ln(1/u)`, `u=c_f(β)`. The remainder — surfaces with holes, wider cross-sections, multi-tube splittings — carries the **signed** products `Wg_γ` of 6j/Weingarten symbols. That signed tail is what decides whether `w₁` survives.

---
## §2 — The candidate inequality: signed-recoupling dominance (precise, with hypotheses)

**Definition (recoupling tail).** Normalize the connected surfaces by the leading tube; each surface `γ` contributes `w(γ)= (Π character coeffs)·(Π 6j/Wg_γ)`. Set
```
𝒯(u) := Σ_{γ≠min, connected} |Wg_γ/Wg_min| · u^{|γ|−|min|}     (each |6j|,|Wg-normalized| ≤ 1).
```
`𝒯(u)` is the *absolute* signed remainder as a fraction of the residue; `|γ|−|min|` is the extra plaquette cost.

**Inequality (Rep-Fourier dominance).** For every `β`, every `T⁴_L`, and all strip lengths `t`:
```
        |A_f|²  =  f₁²  ≥  c₁(β) u^p · (1 − 𝒯(u))²₊ .                        (★)
```
*Proof.* `A_f=A_f^{(0)}(1+Δ)`, `|Δ|≤𝒯(u)` (triangle inequality on the residue-normalized signed tail; `|6j|≤1` bounds each recoupling). `|1+Δ|²≥(1−|Δ|)²`. Divide out `e^{−m_f t}`. ∎

**Hypotheses.** (H1) `u=c_f(β)<1` (small activity). (H_coh) `𝒯(u)≤1−ρ`. Under both, `f₁²≥ρ²c₁u^p=:r₀>0`, and — because the tube weight is *local* and the tail is graded by *cost* — the bound is **`L`-uniform** (KP: the criterion is per-plaquette, identical on every `T⁴_L`).

**§2.2 Why a *lower* bound needs isolation — the Bochner/Hankel form of (★).** Re-expand `C(t)` at integer times as a moment sequence, `a_n=C(nτ)=∫x^n dμ`, `x=e^{−sτ}∈(0,1]`, `μ≥0` (Bochner–Hausdorff: complete monotonicity ⇒ Hankel PSD, total positivity). The atom is a mass at the largest point `x₁=e^{−m₁²τ}`. Any finite-moment, Christoffel/Prony lower bound on `μ({x₁})` has the forced shape
```
   μ({x₁})  ≥  [ Christoffel/Prony functional of a₀…a_{2N−1} ]  −  μ((x₁−ε, x₁)) ,
```
and the second term (*near-edge mass*, what a would-be continuum just under the atom contributes) must be bounded **above** to leave a positive floor. Bounding it is *exactly* "no mass between `m₁` and the next level" = isolation `I2`. Equivalently, `ν` is recovered by the **Fredholm determinant** `D_O(z)=det(1−zT_O)`, `z₁=e^{−am₁²}`; `f₁²=−(1/z₁)N_O(z₁)/D_O′(z₁)` — the same tool the corpus recorded as "看得见原子" (`INVENT-NEW-IDENTITY`, A3). **So the "Bochner lower bound" and the "isolation gap" are one statement in two coordinate systems** (`=GAP-G2-edge` residue–level-gap, re-derived in Fourier coordinates). This is the precise reason the species cannot manufacture a floor from positivity alone: the positivity is `μ≥0` (upper bounds), and the floor is the isolation residual.

**§2.3 What is *new* over prior rep-theory probes.** `H1-T9-spin-network` concluded "spin-network = character expansion; selection rules give `f₁²=0`, never `≥r₀`." That "never" is here made **quantitative**: it is exactly `𝒯(u)≥1−ρ` failing, i.e. the signed tail cancelling. `G2C-MINWEIGHT-RIGIDITY` showed recoupling gives the *sign* not the *floor*. (★) shows the sign is `f⊗f̄→1=+1` and the floor is precisely the dominance of that `+1` over the signed 6j tail — a concrete, computable criterion, not a null statement.

---
## §3 — Three-regime test (honest)

| Regime | `u=c_f(β)` | Signed tail `𝒯(u)` | Does (★) give `>0`? | Verdict |
|---|---|---|---|---|
| **Strong coupling** `β<1/48` | `≈β≪1` (SU(2): `u=0.010`@β=.02, `0.025`@β=.05) | `𝒯(u)≤K u`, `K≤CN` (KP/penrose) `<1` | **Yes** — `f₁²≥ρ²c₁u^p>0`, `V,L`-uniform | **Reproduces P3** (consistency ✓). The grading makes the single tube the unique leading term; `f⊗f̄→1` positive. |
| **Free field** (free Maxwell, `f₁²=0`) | character exp. non-convergent / `g→0` | — | **Yes → 0** | The minimal "tube" is **not a bound state**: `trF²=F·F` is a two-gluon operator, `⟨single | trF² | 0⟩=0` ⟹ `A₁=0`. Numerically the free 2-particle **cut** edge-weight `C(t)e^{s₀t}=½√π·t^{−3/2}→0`, while a genuine atom's stays `=w>0`. **Acid test PASSED** (`w7/verify_signedfourier.py`). |
| **Crossover** `β≈2.5` (Savvidy) | `u=c_f(2.5)≈2.04>1` (SU(2) measured; `c_f>c_0`!) | `𝒯(u)=O(1)`, **uncontrolled** | **No** | `u>1` ⇒ **outside the small-activity radius** `uK<1`; the character expansion no longer converges, `𝒯(u)<1` is a `O(1)`-cancellation question. Restoring it ⟺ the minimal fundamental tube is a **non-cancelling isolated bound state** = `I2` = `m_g>0` = **H1**. **First (and only) break: here.** |

**Reading the crossover column.** The corpus repeatedly found 6j/Weingarten sign-indefiniteness on `SU(N≥3)` (`BS-KERNEL-SIGN`, `T5-signed-weight`, `PROOF-QUANTUM-GROUP-DEFORM`): the 6j symbols have **no termwise sign**, so `𝒯(u)<1` has no representation-theoretic proof — only (a) the small-activity absolute bound (fails at `2.5`), or (b) a *dynamical* input that the tail doesn't cancel: the fundamental tube must be a stable, isolated pole below the two-tube threshold, i.e. **`m_g>0`**. This is the same wall as `BS-RP-STRUCTURE`/`BS-THRESHOLD`: the recoupling is the color-singlet projection, but "continuous → continuous" (composite of a continuum is not an atom), so an isolated atom requires **separate attraction**, and the threshold that attraction must beat is `2m_g`. `𝒯(u)<1` at `2.5` ⟺ that threshold is lifted ⟺ the gap.

**Where (★) *first* breaks (single regime):** going up in `β` from the proven region, it is the **crossover `β≈2.5`** — specifically at the point `u=c_f(β)` crosses `1` (SU(2): between `β=1.0`, `u=0.54`, and `β=2.5`, `u=2.04`). Below that, `u<1` and the absolute tail bound can still work; above, it is a cancellation problem = H1.

---
## §4 — Five-point scorecard

| # | Condition | Verdict | Why |
|---|---|---|---|
| 1 | **Nonperturbative** | **PASS** | Character/spin-foam expansion, RP measure; no loops, no `1/N` in the statement. |
| 2 | **Sign-controlling LOWER bound `m_eff²(k→0)≥λ`** | **FAIL** (direction correct) | (★) is a genuine *lower* bound on the residue — the rare, correct direction — but it does **not** bound `m_eff²`; converting needs `χ≤1/λ` (G2C), so no `λ` is produced across the window. |
| 3 | **Non-circular** | **FAIL** | (H_coh) `𝒯(u)<1` at `2.5` ⟺ `I2` isolation ⟺ `m_g>0` ⟺ H1; (★)'s floor is assumed exactly where it is needed. |
| 4 | **Non-SUSY + fails free field** | **PASS (letter) / hollow on strict read** | The tool is silent (gives `0`/inapplicable) on the free field: `u→1` (`m_f→0`), and the bound-state residue `c₁→0`. Numerically the free cut edge-weight `→t^{−3/2}`. **But** `c₁>0` *is* the sought atom, so the discrimination is **imported, not produced** — this is condition 3's circularity, not an independent pass. |
| 5 | **`(L,a)`-uniform** | **FAIL** | `L`-uniform *at strong coupling* (KP per-plaquette criterion); fails `a→0` and `β→2.5` (`u>1` leaves the radius). |

**Score: 2/5** on the corpus scoring convention (direction-correct lower bound, free-field-silent, but circular — the standard "2/5" shared with BS-RP-STRUCTURE, SUPPORT-NOT-ATOM, DARBOUX, TRANSMUTE), **1/5 on the strict red-team read** where condition 4's free-field silence is judged *imported* (via `c₁>0` = the atom) rather than earned. The two independent adversarial probes `SIGNED-FOURIER-CORE` / `SIGNED-FOURIER-BLIND` take the strict read (1/5, "dead"); the honest spread is **1–2/5**, and the disagreement is *only* whether "correctly silent on the free field because the discriminating input is the atom itself" is a pass or already the circularity. No `5/5`.

The bound is real, direction-correct, free-field-safe, and reproduces the proven theorem — but at the one physical point it needs, `𝒯(u)<1` is the mass gap.

---
## §5 — Sharpest line

The representation-ring / signed-Fourier species does genuine, precise work and then hands back the wall. The Källén–Lehmann measure **is** a non-negative measure on the representation ring (`w_r=|A_r|²≥0`, Bochner/Godement = reflection positivity); the atom is the lightest coefficient; and the only thing standing between `w₁≥0` and `w₁≥r₀>0` is whether the **minimal fundamental tube beats the signed 6j/Weingarten recoupling tail** — the inequality `f₁²≥c₁u^p(1−𝒯(u))²₊`. That inequality is *provably satisfied* at `β<1/48` (where `𝒯=O(u)`, P3) and *provably vacuous* on the free field (where `A₁=0`), exactly as a valid tool must behave. At `β≈2.5` it is undefined: `u=c_f≈2>1`, the small-activity radius that converts the signed tail into an absolute bound is gone, and `𝒯(u)<1` is no longer analysis but dynamics — *the fundamental tube is an isolated, non-cancelling bound state*, i.e. `m_g>0`. The species certifies the atom as a **coherent, sign-definite closure `f⊗f̄→1` whose magnitude is floored by the recoupling network not cancelling** — coherence, not disorder, exactly as advertised — and that non-cancellation is, at the physical point, the mass gap in representation-theoretic clothing.

---
## §6 — Exact stuck step

```
Need:   𝒯(u) ≤ 1 − ρ   (ρ>0)   for β ∈ [β_SC, 2.5],  uniformly in (L, a)
        𝒯(u)=Σ_{γ≠min}|Wg_γ/Wg_min|u^{|γ|−|min|}   =   signed 6j/Weingarten recoupling tail.
Proven β<1/48:  𝒯(u) ≤ K·u < 1  (KP absolute convergence; u≈β≪1)  ⟹ f₁²≥r₀>0  [= P3].
Fails  β≈2.5 :  u = c_f(β) ≈ 2.04 > 1  (SU(2) measured, c_f>c_0)  ⟹ uK≫1; KP dead.
        The remaining content of 𝒯(u)<1 is a O(1) cancellation among sign-indefinite 6j's;
        NO representation-theoretic sign survives (SU(N≥3) 6j is sign-indefinite, no termwise sign).
Equivalent datum:  the minimal fundamental tube is an isolated, non-cancelling bound state
        below the two-tube (string-breaking) threshold:  Λ − m₁² > 0,  Λ ≈ 2 m_g.
        ⟺  isolation I2   ⟺  constituent mass gap m_g > 0   ⟺  H1   ⟺  the mass gap.
```
The missing input is **not** a representation-theoretic lemma (the ring is fully known: positive LR fusion, signed 6j); it is the **dynamical isolation / `m_g>0`** that fixes the *magnitude* of the recoupling non-cancellation. Positivity/recoupling fix the **sign**; only confinement fixes the **floor**.

---
## §7 — Numeric acid-test appendix (`w7/verify_signedfourier.py`, `py -3.12`, numpy/scipy)
- SU(2) `c_j(β)>0` for `β∈{0.02,0.1,0.5,1.0,2.5}`, `c_0=2I₁(β)/β` exact; fundamental `j=½` minimal Casimir.
- **Crossover signature (measured):** at `β=2.5`, `c_{½}=2.04 > c_0=2.01` — `u>1`, the KP radius is long exceeded. (The proven radius is `β<β₁=O(1/N)`, `1/48` in SZZ; `u≍β` there, so `uK<1` holds only at `β≲O(1/N)`.)
- **Atom vs cut (acid test):** edge-weight `C(t)e^{s₀t}` = const `=w>0` for an atom; `=½√π·t^{−3/2}→0` for the free 2-particle cut. `f₁²=0` on the free field, correctly.

---
*Refs.* Drouffe–Itzykson, Phys. Rep. **38C** (1978) 133 (positive character coefficients); Itzykson–Zuber, J. Math. Phys. **21** (1980) 411; Godement/Bochner positive-definite functions on groups; Osterwalder–Seiler, Ann. Phys. **110** (1978) 440; Kotecký–Preiss, CMP **103** (1986) 19; Shen–Zhu–Zhu, CMP **400** (2023) 805, arXiv:2204.12737; Collins–Śniady, CMP **264** (2006) 773 (signed Weingarten); Varshalovich–Moskalov–Khersonskii, *Quantum Theory of Angular Momentum* (6j signs); Caron-Huot–Pokraka–Zahraee, JHEP **01** (2024) 195 (positivity alone ⇒ residue only `≥0`). Internal: PROOF-P3-UNCONDITIONAL (`f₁²≥c₁u^p(1−O(u))`, `f⊗f̄→1`), PROOF-P1P2, PROOF-ISOLATION (I2 open), PROOF-F1SQ (free `:φ²:` `f₁²=0`), GAP-G2-edge (residue–level-gap lemma; moments ⇒ upper), INEQ-R10-sumrule-positive (only upper sum rules; `arXiv:2309.04472`), G2C-MINWEIGHT-RIGIDITY (recoupling = sign not floor; GFF), H1-T9-spin-network (character expansion; no rep-theoretic floor — made quantitative here), PROOF-HEAT-KERNEL (character/heat-kernel gap = fiber, not 0⁺⁺), PROOF-QUANTUM-GROUP-DEFORM (root-of-unity positivity is truncation-bound; Weingarten sign returns), T5-signed-weight, BS-KERNEL-SIGN / BS-RP-STRUCTURE / BS-THRESHOLD (recoupling sign-indefinite; composite→composite; threshold `2m_g`), INVENT-NEW-IDENTITY (Fredholm `D_O(z)` sees the atom = same object), SUPPORT-NOT-ATOM, PROOF-SKELETON (G2.5/G2.6); sibling red-teams `SIGNED-FOURIER-CORE`, `SIGNED-FOURIER-BLIND` (independent 1/5 "dead" — strict read); `w7/verify_signedfourier.py`.
