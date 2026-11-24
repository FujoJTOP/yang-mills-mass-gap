# THIMBLE-BUILD — the Picard–Lefschetz / resurgence species, built as a precise inequality: the atom `f₁²` is the **modulus square of a coherent thimble sum**, `f₁²=|\sum_\sigma n_\sigma a_\sigma e^{i\,\mathrm{Im}\,S_\sigma}|² ≥ \cos\delta\,(\sum_\sigma a_\sigma)²`; the phase is *unobservable in a single thimble* but **survives the modulus square as the inter-saddle interference weight `\cos\delta`**, and its control at `β≈2.5` is exactly the Stokes constant `=⟨\mathrm{tr}F²⟩` `=H1`

> 2026-11-23 · 4D pure SU(N) YM · Constructive forge of the **Picard–Lefschetz thimble / resurgence** species for the atom `f₁²=|⟨0|\mathrm{tr}F²|G⟩|²`, *built as an inequality with stated hypotheses* — not prose.
> **判词 (verdict).** The construction is real, precise, and **materially better-directed than the dead 6th species**: `f₁²` is a genuine **modulus square**, hence a **coherence (Schur-convex) functional** — *maximized* by phase alignment, *destroyed* by phase spreading — which is exactly the `total ≤ atom` direction that `SIGNED-FOURIER-CORE` proved the representation-ring species could **never** supply (its face was Schur-concave, `atom ≤ total`). The one-line result is `f₁² ≥ \cos\delta\cdot R²`, `R=\sum_\sigma a_\sigma`. It **(i) reduces to P3 at `β<1/48`** (`\delta=0`, real contour, `R=a_{\rm tube}`, consistency); **(ii) correctly returns `0` on the free field** (single gapless saddle: the isolated-mode residue `a_0=0`, no atom — acid test PASSED, but the pass *is* the gaplessness); **(iii) breaks at the crossover `β≈2.5`**, not on direction but on the **phase**: the Savvidy saddle carries `\mathrm{Im}\,S_\sigma=c(gB)²V>0`, the bion/instanton tower spreads `\mathrm{Im}\,S` densely on the circle, so `\cos\delta` has **no positive lower bound** and the bound degrades to the free-field value `f₁²≥0`. The phase `e^{i\,\mathrm{Im}\,S}` dies in `|A_1|^2` only as an *overall* phase (single thimble: unobservable); it **survives as the relative phase**, and requiring `\cos\delta≥\cos\delta_0>0` is *literally* a one-sided bound on the Stokes constant `σ_{\rm Stokes}=⟨\mathrm{tr}F²⟩=Λ=H1`. **Score 2/5** (charitable; 1.5/5 strict, §5). Sharpest: the thimble modulus-square supplies the missing **direction** (`total` bounded by `atom`) that the 6th species lacked — and charges for it in the **Stokes/condensate** that the 6th species also could not pay.

---

## §0 · The one-line inequality

```
Complexified action S(U), links U∈SL(N,ℂ), saddles {φ_σ}, thimbles J_σ, intersection numbers n_σ.
Isolated 0⁺⁺-mode residue of saddle σ:   a_σ e^{iφ_σ},   a_σ ≥ 0,   φ_σ := Im S_σ + arg(det' M_σ) + arg(overlap).

  f₁²  =  | A₁ |²  =  | Σ_σ n_σ a_σ e^{iφ_σ} |²
       =  Σ_σ a_σ²  +  2 Σ_{σ<σ'} a_σ a_{σ'} cos(φ_σ−φ_{σ'})
       ≥  (cos δ) · ( Σ_σ |n_σ| a_σ )² ,      cos δ := min_{σ,σ'∈supp} cos(φ_σ−φ_{σ'}) .        (★)
```
**Hypotheses.** (A) the thimble decomposition converges to the physical reflection-positive correlator, uniformly in `t`, so residue extraction commutes with `Σ_σ`; (B) each contributing saddle carries the isolated `0⁺⁺` mode with real residue modulus `a_σ≥0`; (C) **phase coherence**: `cos δ ≥ cos δ₀ > 0`; (D) **residue floor**: `R:=Σ_σ|n_σ|a_σ ≥ R₀>0`. Under (A)–(D): `f₁² ≥ cos δ₀ · R₀² > 0`. **Hypotheses (B) at the free field and (C)(D) at `β≈2.5` are the whole game** (§4).

---

## §1 · The thimble construction: saddles, residue, and which combination carries `f₁²`

**1.1 Complexification and the saddle set.** On `T⁴_L` (spacing `a`) with Wilson action `S_W=−(β/N)Σ_p Re tr U_p`, continue the links `U_e∈SU(N)` to `U_e∈SL(N,ℂ)` and `S_W(U)=−(β/2N)Σ_p(\mathrm{tr}\,U_p+\mathrm{tr}\,U_p^{-1})` holomorphically (Cristoforetti–Di Renzo–Mukherjee, Witten). Critical points `∂S_W/∂U=0` are, on the lattice, the higher-coupling analogues of the continuum **Savvidy saddles**: constant abelian `F_{12}=B` embedded in a Cartan direction, plus the trivial `φ_0` (`U=1`), plus the **multi-bion / multi-instanton tower** `φ_{n}` obtained by stacking. The continuum Hessian on transverse 1-forms is exact (`PROOF-SAVVIDY-EXACT-SIGN`, `PROOF-SAVVIDY-REFORMULATION-CONVEX`):
```
H  =  −D² + 2\,ad F_B ,        spectrum on constant B:   λ_{n,s}=k_∥² + (2n+1)g|B| + 2s\,gB ,
λ_min = −g|B|  at  n=0, s=−, k_∥=0   (softest at k_∥=0),
N_neg/V = c_N (gB)²/8π²  (c_{SU2}=1, c_{SU3}=3),   Tr(H_−) ≤ 8C_N ∫|F|³ .
```
This is the species' saddle data verbatim: **negative transverse-magnetic modes, depth `−g|B|`, count `∝(gB)²V`, all monotone in `|B|`, softest at `k_∥=0`**. At `β=2.5` the *measured* lattice Hessian reaches `min=−1.12β` with `22–49%` negative directions (`W6-BE-at-2.5`), i.e. the Savvidy instability is present and physical (gauge-quotient survived, `OTHER-SAVVIDY-SIGN`).

**1.2 The thimble decomposition.** With the flow `\dot U=\overline{∂\bar S}` descending from each `φ_σ`, the stable thimble `J_σ` has `\mathrm{Im}\,S≡\mathrm{Im}\,S_σ` (the flow conserves it) and `|e^{-S}|=e^{-\mathrm{Re}\,S}>0`. For the `\mathrm{tr}F²` correlator `C(t)=⟨O(t)O(0)⟩_c`,
```
C(t) = Σ_σ n_σ e^{-i\,\mathrm{Im}\,S_σ}\, e^{-\mathrm{Re}\,S_σ}
             ∫_{J_σ} e^{iω_σ}\, O[U(t)]O[U(0)]\, e^{-(S-S_σ)}\,dμ_σ
     = Σ_σ n_σ e^{i\,\mathrm{Im}\,S_σ}\, e^{-\mathrm{Re}\,S_σ}\, c_σ(t)\, [\,1+O(\hbar)\,] ,
```
with `ω_σ=\arg\det J_σ` the residual Jacobian phase. The Källén–Lehmann measure `ν≥0` (`C(t)=∫e^{-mt}dν`, `f₁²=ν({m₁²})`) is recovered by residue extraction. **The atom is the coefficient of the isolated exponential `e^{-m_1 t}`.**

**1.3 Which combination carries `f₁²`.** Taking the residue at the isolated `0⁺⁺` pole `m₁` of the *σ-thimble term* gives `c_σ(t) ⊃ a_σ e^{i\,\mathrm{Im}\,S_σ}e^{-m_1 t}` where `a_σ=e^{-\mathrm{Re}\,S_σ}|\det' M_σ|^{-1/2}|\langle G|O|\phi_σ\rangle| \ge 0` (the bracket absorbs the Jacobian and overlap phases into `φ_σ`). Hence
```
A₁ = lim_{t→∞} e^{m_1 t} C(t) = Σ_σ n_σ a_σ e^{iφ_σ} ,      φ_σ = Im S_σ + (Jacobian+overlap phase) .     (♦)
```
**`f₁²=|A₁|²` carries the modulus square of a coherent thimble sum.** *Not* `Z` (no modulus square), *not* a single determinant (no coherence), *not* a positive-type functional (Schur-concave). It is the one object in the campaign built to detect **coherent constructive interference of the saddle tower at the isolated pole** — precisely the "atom is a coherence, not a sum" functional.

---

## §2 · The candidate inequality, with proof

**Theorem (conditional thimble lower bound — the algebra is unconditional; the hypotheses are the physics).** *Assume (A)(B) of §0. Let `\mathrm{supp}` be the saddles of non-negligible weight in (♦), `a_σ≥0`, `φ_σ∈ℝ/2π`. Then*
```
f₁²  =  Σ_σ a_σ²  +  2 Σ_{σ<σ'} a_σ a_{σ'} cos(φ_σ−φ_{σ'})  ≥  (cos δ)(Σ_σ a_σ)² ,   δ := max_{σ,σ'} |φ_σ−φ_{σ'}| mod 2π ,     (★)
```
*with the convention that if `cos δ<0` the RHS is `≤0` and the bound is trivial (`f₁²≥0`).*

*Proof.* `f₁²=|Σ a_σ e^{iφ_σ}|²=Σ_{σ,σ'}a_σ a_{σ'}cos(φ_σ−φ_{σ'})`. By hypothesis `cos(φ_σ−φ_{σ'})≥\cos\delta ≥0` for all pairs. Then
`f₁² ≥ Σ a_σ² + 2\cos\delta\,Σ_{σ<σ'}a_σ a_{σ'} = \cos\delta\,(Σ a_σ)² + (1-\cos\delta)Σ a_σ² \ge \cos\delta\,(Σ a_σ)²`,
using `1-\cos\delta≥0` and `Σ a_σ²≥0`. ∎ (Verified numerically, 0 violations in 46 676 random towers with `\cos\delta>0`; §8.)

**Deficit form (equivalent, shows what eats the atom).**
```
f₁²  =  (Σ_σ a_σ)²  −  2 Σ_{σ<σ'} a_σ a_{σ'} (1 − cos(φ_σ−φ_{σ'}))  .        (★')
```
The **absolute maximum** `(Σa_σ)²` is attained iff all phases align (`φ_σ≡\phi`, single global phase); every unit of phase spread removes `2a_σ a_{σ'}(1-\cos\psi)` from the atom. **The atom is a coherence functional**: it is Schur-**convex** in the phase distribution (maximised by concentrating all phase, minimised by spreading) — the exact opposite of a positive-type/Bochner functional.

**Sharpening (single real saddle ⇒ no phase problem).** If the dominant saddle is real and unique up to gauge (`φ_σ=0`, `\mathrm{Im}\,S=0`), then `f₁²=a_σ²=(Σa_σ)²`, bound `(★)` is an equality. The phase problem is **intrinsically a multi-saddle / near-Stokes phenomenon**, invisible to single-thimble dominance.

---

## §3 · The crux: does the phase `e^{i\,\mathrm{Im}\,S_σ}` survive the modulus square `f₁²=|A₁|²`?

This is exactly where the 6th species (`SIGNED-FOURIER-CORE`) died, so answer it case by case.

**3.1 The absolute phase dies — unobservable.** `f₁²=|A_1|²` kills a *common* phase: if `A_1=e^{i\phi}\sum a_σ` (all `φ_σ` equal), then `f₁²=(Σa_σ)²`, phase gone. So `e^{i\,\mathrm{Im}\,S_σ}` does **not** appear in `f₁²` as an overall factor. *A single thimble contributes its residue squared, full stop* — the free-field/strong-coupling regimes are phase-free and the bound is tight (§4).

**3.2 The relative phase survives — this is the lever.** With ≥2 saddles, `f₁²` keeps `\cos(φ_σ−φ_{σ'})`. The `\mathrm{Im}\,S` differences are **not** gauge-removable (they are the trans-series exponent between distinct complex saddles). They enter `f₁²` as genuine **interference weights** `\cos\delta∈[−1,1]`:
```
 phase alignment  δ=0     ⇒  f₁² = R²          (maximal, constructive)
 phase spread     δ=π/2   ⇒  f₁² ≥ 0           (the bound collapses to the free-field value)
 phase spread     δ>π/2   ⇒  bound is trivial  (cosδ<0; f₁²≥0 only)
```
So `\mathrm{Im}\,S` survives the modulus square as **degradation**, never as an overall sign. The 6th species' problem was different and opposite: its "signed" content (Weingarten alternation) lived in a *real* surface sum (no oscillation), and its usable face (`c_R≥0`) was Schur-concave ⇒ could only give `atom ≤ total`, the **wrong direction**. Here the direction is correct (Schur-convex `total ≤ atom`), and the price is a **genuinely oscillatory** phase `\mathrm{Im}\,S_σ ≠0` at the Savvidy saddle.

**3.3 Where the Stokes / condensate ambiguity enters — a one-saddle-at-a-time verdict.** The thimble sum is a **trans-series**: asymptotically, its lateral sums on the two sides of a Stokes line differ by
```
ΔA₁  =  σ_Stokes · e^{−S_sav} · A₁^{(J_σ*)} + … ,        S_sav = 8π²/g² ,
```
with `σ_Stokes` undetermined on `ℝ⁴` (`RESURGENCE` §2; Maiezza–Vasquez 2405.01639). Crossing the Stokes line **shifts `φ_σ*` by the Stokes phase** — i.e. it moves `\cos\delta`. Therefore
```
   Hypothesis (C) "cos δ ≥ cos δ₀ > 0  uniformly"  ⟺  a one-sided bound on the Stokes phase
                                                    ⟺  a one-sided bound on σ_Stokes
                                                    ⟺  a one-sided bound on ⟨tr F²⟩ = Λ = H1 .
```
Concretely, the Savvidy saddle's one-loop effective action has `\mathrm{Im}\,V_eff ∝ (gB)²>0` (Nielsen–Olesen negative mode) and the bion tower `φ_n` accumulates phase `n·ψ_bion`; as `n→∞` or `V→∞` the phases `φ_σ` become **dense on the circle**, so `min_{σ,σ'}\cos(\phi_σ-\phi_{σ'}) = −1` and `\cos\delta` has **no positive lower bound**. Restoring `\cos\delta≥\cos\delta_0>0` *is* fixing the condensate. **This is the precise statement of the undetermined term** — and per the task's charge, it **is** the valuable output: the atom is a coherent modulus square whose coherence, at the physical point, is the Stokes constant.

**3.4 Bonus: the failure is `V`-extensive, so it also fails (L,a)-uniformity.** For the homogeneous Savvidy saddle `\mathrm{Im}\,S_σ \propto (gB)²V` is **extensive**; mod `2π` it wraps as `V→∞`. Even if `\cos\delta>0` at small volume, there is **no `L`-uniform floor** — the coherence hypothesis is violated precisely in the thermodynamic limit where the gap must be certified. (This is the concrete mechanism behind condition 5's failure, §5.)

---

## §4 · Three-regime test (honest)

| Regime | Saddle content | `a_σ`, phases | Does (★) give `>0`? | Verdict |
|---|---|---|---|---|
| **Strong coupling `β<1/48`** | Physical contour **real**, `\mathrm{Im}\,S≡0`; the "thimble sum" is the positive polymer sum, leading term = single fundamental tube (`PROOF-P3`) | `a_{\rm tube}=a_0u^{p}(u)^{t/a}`, `u=c_f(β)`, `φ=0` (all `φ_σ=0`) | **Yes** — `\delta=0`, `f₁² ≥ R² = a_{\rm tube}² = r₀>0`, `V,L`-uniform (`KP`) | **Reduces to P3** (consistency ✓). Phase machinery idle (`δ=0`); like `RESURGENCE` §1, adds nothing here. |
| **Free field `g=0`** (`f₁²=0`) | Single saddle `φ_0` (`A=0`); `H=−D²⪰0`, **gapless**; `ν` = two-photon **cut** | `a_0` := **isolated-mode** residue = `0` (no isolated pole to take residue of; `\mathrm{tr}F²=F\!\cdot\!F` is a 2-gluon operator, `⟨\text{single}|O|0⟩=0`) | **Returns 0** — `R=\sum a_σ=0` ⟹ `f₁²≥0`, in fact `=0` | **Acid test PASSED** — *and the pass is the gaplessness*: no isolated mode ⇒ `a_0=0`. **Vacuity trap bypassed only by the isolated-mode definition** (§4.1). |
| **Crossover `β≈2.5`** (Savvidy) | Trivial `φ_0` (`\mathrm{Im}\,S=0`) + Savvidy `φ_B` (`\mathrm{Im}\,S∝(gB)²V`) + bion tower `φ_n` (`\mathrm{Im}\,S∝n`) | `a_σ>0` (semiclassical `e^{-S_sav}`), but **phases spread densely** on `ℝ/2π` | **No** — `\cos\delta` has no positive lower bound (`→−1` as tower/`V` grow); (★) degrades to `f₁²≥0` | **First (and only) break: the phase `\cos\delta`.** Restoring it = fixing `σ_Stokes=⟨\mathrm{tr}F²⟩=H1` (§3.3). |

**4.1 The free-field pass and its trap (honesty note).** The pass hinges on `a_0` being the **isolated-mode** residue. If one instead reads `a_0` off the one-loop determinant of the Gaussian saddle (treating the whole continuum as "the mode"), one gets `a_0>0` and (★) **falsely** returns `f₁²>0` on the free field — **VACUOUS**, the same trap that killed the 6th species' benign face. The correct object is the residue of an *isolated* pole; the free field has none. So the acid test is passed, but — exactly as `RESURGENCE` §3 — **the pass is a restatement of the failure to manufacture a floor from a gapless saddle**. The free-field `0` is not produced by thimble dynamics; it is produced by the absence of an isolated mode. This is why §5 scores (4) a pass but (2) only a half.

**4.2 Reading the crossover column.** The Savvidy background is *not* the only contributor, which is the whole difficulty: `φ_0` (real, `\mathrm{Im}\,S=0`) and `φ_B` (complex, `\mathrm{Im}\,S∝(gB)²V>0`) and the tower `φ_n` carry **pairwise relative phases** `ψ_{0,B}=\mathrm{Im}\,S_B`, `ψ_{n,m}∝(n−m)`. As `B` ranges over the flux lattice and `n` over the tower, `{ψ}` is dense mod `2π`, so `\cos\delta:=\min\cos\psi\to−1`. There is **no thimble selection rule** — no positivity of `\cos\delta` — because the phases are the trans-series exponents, and their sum is the Borel/Stokes-ambiguous series (`RESURGENCE`). The bound's collapse `(★)→(f₁²≥0)` is not a loose estimate; it is the *generic* value once phases decohere.

**Where (★) *first* breaks (single step):** going up in `β`, it is the point where the Savvidy saddle's **`\mathrm{Im}\,S` first becomes `O(1)` mod `2π`** — the **Stokes line** of the `e^{-S_sav}` trans-series — at the **crossover `β≈2.5`**. Below it phases cluster (`\cos\delta>0` possible); on/above it the bion tower fills the circle and `\cos\delta` goes negative. That step is `σ_Stokes`, the condensate.

---

## §5 · Five-point scorecard

| # | Condition | Verdict | Why |
|---|---|---|---|
| 1 | **Nonperturbative-rigorous** | **PARTIAL (½)** | (★)–(★') are exact algebra on the trans-series coefficients (no `\hbar` truncation error); the *representation* is genuinely nonperturbative (saddles, Stokes). But the equality "thimble sum = physical RP correlator uniformly in `t`" (hyp. A) is **unproven at `β≈2.5`** (`T5`; `PROOF-THIMBLE` "Cauchy invariance" only gives *same theory*, not *controllable residue*). Rigorous at SC and free field. |
| 2 | **Sign-controlling LOWER bound** | **PARTIAL (½)** | (★) is a genuine **lower** bound in the **Schur-convex** direction — the `total ≤ atom` inequality that `SIGNED-FOURIER-CORE` proved representation-positivity **cannot** give. *New and correct direction.* But it is sign-controlling only where `\cos\delta>0` is *provable* (SC: yes, `δ=0`; free: vacuous; **crossover: no**). Direction achieved in the formula, magnitude not controlled at the target. |
| 3 | **Non-circular** | **FAIL (0)** | Hypothesis (C) `\cos\delta≥\cos\delta_0>0` at `β≈2.5` is, by §3.3, a one-sided bound on `σ_Stokes=⟨\mathrm{tr}F²⟩=Λ` — the condensate — i.e. H1. The floor is **assumed exactly where it is needed**. Same circularity as `RESURGENCE` #3. |
| 4 | **Non-SUSY + fails free field** | **PASS (1)** | Pure YM, no SUSY anywhere. Free field: single gapless saddle ⇒ isolated-mode residue `a_0=0` ⇒ `R=0` ⇒ `f₁²≥0`, and `=0`. **Not vacuous** (it returns the true `0`), *provided* the isolated-mode definition (§4.1). The positive-type comparison never appears. |
| 5 | **`(L,a)`-uniform** | **FAIL (0)** | `\cos\delta` has **no positive `L`-uniform floor**: `\mathrm{Im}\,S∝V` is extensive and wraps mod `2π`; the atom is usc and may drop to `0` at the crossing (`PROOF-TRANSPORT-ANALYTIC` K1). No `a→0` control of the Stokes phase. |

**Score: 2/5** (conditions #4 full, #1 and #2 halves; #3, #5 fail). **Honest spread 1.5–2/5.** The *strict* read takes 1.5/5: condition (4)'s free-field pass imports the isolated-mode definition (`a_0=0` is *the atom*, `=H1` at general `β`), so the discrimination is partly **imported, not produced** — the same tension flagged in `SIGNED-FOURIER-BUILD`'s 1-vs-2 split. The charitable read gives 2/5 because, unlike the 6th species, the thimble formula **does supply the missing Schur-convex direction** and **is exactly tight in the two anchored regimes** (`δ=0` at SC ⇒ P3; `R=0` at free ⇒ 0). Either way: **not 5/5**, and the gap is the **Stokes phase**, not the direction.

*Comparative (6th species vs thimble).* `SIGNED-FOURIER-CORE` died on the **direction**: Bochner/positive-type is Schur-**concave** (`atom ≤ total`), so it can never give `f₁²≥r₀`. The thimble modulus square is the first species to put the **direction right** (`f₁²=|\sum|²` is Schur-**convex**, `total ≤ atom`). What it charges instead is the **phase**: the Schur-convex coercivity is coherent-interference-limited, and the coherence `\cos\delta` is the Stokes constant. So the thimble species and the 6th species fail on **two complementary walls of the same H1** — direction vs magnitude/phase. Neither crosses.

---

## §6 · Sharpest line

*The atom is a modulus square, so it is the campaign's first genuinely **Schur-convex** functional of the saddle tower — `f₁²=|\sum_\sigma a_\sigma e^{i\,\mathrm{Im}\,S_\sigma}|²`, maximised by phase alignment, killed by phase spread. This is exactly the `total ≤ atom` lever that the representation-ring species (`SIGNED-FOURIER-CORE`) proved positivity could never build; on the **direction** it is strictly better. But the modulus square does not make the phase vanish — it makes it **relative**: `f₁²=(\sum a_\sigma)²−2\sum_{σ<σ'}a_\sigma a_{\sigma'}(1-\cos\psi)`, and the interference factor `\cos\delta` **is** the trans-series Stokes phase. A single thimble is phase-free and tight (`δ=0`, reproduces P3; `R=0`, correctly zero on the free field); the crossover is where the Savvidy `\mathrm{Im}\,S∝(gB)²V` and the bion tower `∝n` make `\{\psi\}` dense on the circle, `\cos\delta\to−1`, and the bound collapses to `f₁²≥0`. Requiring `\cos\delta≥\cos\delta_0>0` is a one-sided bound on the Stokes constant — `σ_Stokes=⟨\mathrm{tr}F²⟩=Λ=H1`. The thimble species buys the missing **direction** with the **condensate**; the 6th species lacked both; both are the same gap.*

---

## §7 · Exact stuck step

```
Need (the one missing lemma, [T-COH]):  phase coherence of the isolated-mode thimble tower,
        cos δ  :=  min_{σ,σ'∈supp} cos( Im S_σ − Im S_σ' + Δ(Jacobian+overlap) )
               ≥  cos δ₀ > 0        for β ∈ [β_SC, 2.5],  uniformly in (L, a),
   so that   f₁² = |Σ_σ n_σ a_σ e^{i Im S_σ}|²  ≥  cos δ₀ · R₀²  > 0.

Proven:  SC β<1/48 — Im S ≡ 0 on the real contour ⇒ δ=0, R=a_tube ⇒ f₁² ≥ r₀>0   [= P3].
         Free g=0  — single gapless saddle ⇒ a_0=0 (isolated-mode residue) ⇒ R=0 ⇒ f₁² ≥ 0, =0.
Fails:   β≈2.5   — Savvidy Im S_σ ∝ (gB)²V > 0 (Nielsen–Olesen negative mode) and bion tower Im S ∝ n
                   ⇒ {ψ} dense mod 2π ⇒ cos δ → −1 (no positive floor); (★) → f₁²≥0 (vacuous).
         Keeping δ bounded is, identically,   σ_Stokes ≤ σ₀   ⇔   ⟨tr F²⟩ bounded one side
                   ⇔   the gluon condensate   ⇔   Λ   ⇔   H1   ⇔   the mass gap.
```
The missing input is **not** thimble selection / saddle dominance / a positivity identity (all checked, `PROOF-THIMBLE`), and **not** the direction (the modulus square supplies it). It is **phase coherence of the isolated-mode tower** — a one-sided bound on the Stokes constant. The formal trans-series is
```
f₁²  =  | Σ_σ n_σ e^{−S_σ} (det' M_σ)^{−1/2} ⟨G|O|φ_σ⟩ |²    (formal; sads σ over Savvidy+bion tower),
```
**an exact formal modulus square with an undetermined Stokes constant `σ_Stokes=⟨tr F²⟩` — stated precisely, that is the deliverable**: the atom is *coherent*, and its coherence at the physical point *is* the condensate.

---

## §8 · Numeric acid-check (`py -3.12`, numpy 2.4.6; trivial scalar tower, no lattice job)

- **Bound (★) verified:** 300 000 random towers (`2≤K≤6` saddles, `a_σ∈[0,1]`, `φ_σ∈[−π,π]`); restricted to `\cos\delta>0` (46 676 cases) there are **0 violations** of `f₁² ≥ \cos\delta\,(\sum a_\sigma)²`; `max(rhs−f₁²)=0`.
- **Coherence collapse (2-saddle, `a=(1,0.5)`):** `(Σa)²=2.25`; `ψ=0 ⇒ f₁²=2.25` (tight, `\cos\delta·R²=2.25`); `ψ=π/4 ⇒ 1.957` (bound `1.591`); `ψ=π/2 ⇒ 1.25` (bound `0` — collapse); `ψ=π ⇒ 0.25`. Deficit form `(★')` reproduces `f₁²` **exactly** at every `ψ`. The bound degrades to `f₁²≥0` at `\cos\delta=0` (`ψ=π/2`) — the codimension of the crossover failure.
- **Free-field edge-weight (acid anchor, matches `SIGNED-FOURIER`/`PROOF-F1SQ`):** `C(t)e^{s_0t}=½√π\,t^{−3/2}→0` for the two-gluon cut (atom `f₁²=0`), versus `=\text{const}=w>0` for a genuine isolated pole. The thimble residue `a_0` must track the *pole*, hence `a_0=0` on the free field — the load-bearing definition of §4.1.

---

*Refs.* External: M. Cristoforetti, F. Di Renzo, L. Scorzato, arXiv:1205.3996; G. Aarts et al., arXiv:1303.7204; E. Witten, arXiv:1001.2933 (analytic continuation / thimbles); G. V. Dunne, M. Ünsal, arXiv:1210.2423; Argyres–Ünsal, arXiv:1206.1890 (resurgence, trans-series); Aniceto–Başar–Schiappa, Phys. Rept. **809** (2019); Costin, *Asymptotics and Borel Summability* (2009); 't Hooft, 1977 (renormalon = condensate); Beneke, Phys. Rept. **317** (1999) 1; Maiezza–Vásquez, 2405.01639 (`\sigma`-ambiguity on `ℝ⁴`); Savvidy, PLB **71** (1977) 133; Nielsen–Olesen, NPB **144** (1978) 376; Leutwyler, NPB **179** (1981) 129 (tachyonic/imaginary part `∝(gB)²`); Osterwalder–Seiler, Ann. Phys. **110** (1978) 440 (RP); Shen–Zhu–Zhu, arXiv:2204.12737, CMP **400** (2023) 805 (rigorous `β<1/48`).
Internal: `PROOF-SAVVIDY-EXACT-SIGN` (Hessian `−D²+2adF_B`, `λ_min=−g|B|`, `N_neg∝(gB)²V`, `k_∥=0` softest); `PROOF-SAVVIDY-REFORMULATION-CONVEX`; `W6-BE-at-2.5` (measured `−1.12β`, `22–49%` neg); `OTHER-SAVVIDY-SIGN` (gauge-quotient negativity physical); `PROOF-THIMBLE` (thimbles = relocated phase, not positivity; Cauchy invariance; T1/T2 = H1); `RESURGENCE` (convergent `β<β₁`; Stokes `σ` = condensate; `0↦0` free-field pass = impotence); `SIGNED-FOURIER-CORE`/`SIGNED-FOURIER-BUILD` (Schur-**concave** direction wall; the complementary failure); `PROOF-P3-UNCONDITIONAL` (SC `f₁²≥r₀>0`); `PROOF-F1SQ`/`PROOF-HEAT-KERNEL` (free `f₁²=0`, 2-photon cut); `PROOF-TRANSPORT-ANALYTIC` (K1, `f₁²` usc); `SCHUR-CONVEX`/`EDGE-f1²-attack` (master sign obstruction; this species attacks it from the convex side); `W2-transmuted-gap`, `T5-signed-weight`.
