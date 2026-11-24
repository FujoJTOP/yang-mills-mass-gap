# E2-FLIP-DYNAMICAL — attempt to CONSTRUCT the coupling-dependent, non-Schur, non-regress functional Φ that lower-controls E₁ across the Savvidy crossover; verdict: **IRREDUCIBLE — the "Savvidy-as-lever" direction is structurally impossible, and every Φ collapses to H1**

> 2026-11-24 · Direct attack on the precise **E₂-flip condition** of `E2-REGRESS-CIRCULARITY §6`, using the one place the dynamics is nonperturbatively special — the **Savvidy / crossover structure** — explicitly as a suggested **lever** rather than an obstacle (`CLAUDE.md §5`: "the Savvidy/crossover structure is the place to look").
>
> **Labels (`CLAUDE.md §2`):** **[T]** 真定理 / **[C]** 条件定理 / **[NG]** no-go *for a named class* / **[O]** 精确缺失步 / **[J]** 判断 / **[V-here]** fresh computation this session. External facts are **[R]** (standard/as-reported, not re-verified here). **No QED on anything unproven. This document does NOT prove Clay and does NOT prove that no tool can exist (`∀class ≠ ∀possible`, `KEY-EXISTS §5`).**

---

## Verdict up front

> ### **IRREDUCIBLE. No Φ constructed. The Savvidy-Hessian lever is provably the *wrong direction*: the Hessian depth and the true transfer gap are anti-correlated, so the more Savvidy-unstable the theory, the *smaller* its gap — any monotone functional of the action Hessian is a lower bound on nothing (indeed it goes negative where the gap is positive).**
>
> The single-page reason, four clauses (each a distinct, named, *survivable* failure of one natural Φ — not one repeated excuse):
> 1. **Wrong object (UV fluct. Hessian vs IR transfer gap) — proven in a solvable witness [T]+[V-here].** The action Hessian `Hess V = −D²+2adF_B` is the **second variation of the action**; the gap is the **spectral gap of the conditioned measure (transfer-matrix / Langevin generator)**. In an exactly-solvable model these are shown to be **anti-correlated**: as the wells deepen (Savvidy depth `−gv² → −∞`), the true gap collapses **exponentially** (`≈ e^{−S_inst}`, `S_inst ∝ √g`): `gap/|depth| → 0` over ≈5 orders of magnitude. So "use Savvidy as the lever" is **impossible** — Savvidy negativity is precisely what makes a Hessian-functional *decrease* where the gap is forced *positive*.
> 2. **Wrong direction (conditioning, not sign) — [T].** At the unconditioned level the corpus's exact Stein identity `E_μ[Hess V] = E_μ[∇V⊗∇V] ⪰ 0` holds (`H1-HESSIAN-DIRECT`); the toy confirms it numerically to 3–4 digits. But this expected Hessian is **gap-blind** and even *diverges positively* while the gap shrinks (toy: `E[H] ≈ 2g` as `g→∞`, gap `→0`). The gap lives in the **conditioned** curvature `G″ = E[H]−Var+W = S−D` (`=H1`), one conditioning level above the Savvidy Hessian.
> 3. **Not uniform — [NG].** Savvidy data is **extensive/UV**: negative-mode count `N_neg = c_N(gB)²V/8π² ∝ V`, depth `−g|B| ∝ β` (`PROOF-SAVVIDY-EXACT-SIGN`). A quantity `∝ V` or `∝ β` cannot be an `(L,a)`-uniform lower bound on the IR mass.
> 4. **Non-regress but *circular repackaged* — [O].** The only *direction-correct* Savvidy-derived functional — the **tunnelling/instanton (condensation) action** that *produces* the mass — equals the gap but **computing/controlling it uniformly across the crossover IS the nonperturbative transport machine = H1**. So the direction-correct version is not a flip; it is the missing machine under its standard name.
>
> **The precise incompatibility:** the gap, in a marginal scale-free theory that generates its mass by **tunnelling between unstable configurations** (Savvidy/dimensional transmutation), is a **global, non-local, conditioned** functional of the measure, *not* a function of the local second variation of the action. The E₂-flip condition asked for a functional that is *simultaneously* (i) non-Schur, (ii) non-regress, (iii) dynamical (Savvidy), (iv) acid-test-clean, (v) uniform, (vi) non-circular. Clauses (iii)+(i) push one to the Hessian; clauses (v)+(vi) are violated by the Hessian; clause (iv) alone is satisfied but buys nothing. **The four constraints are mutually inconsistent for the Savvidy lever — this is a genuine incompatibility, not (only) a failure of this attempt** (scoped to the named class "monotone/continuous functionals of the action's local second variation").

---

## §0 The flip condition, restated with the exact object

From `E2-REGRESS-CIRCULARITY §6`, a flip needs a functional `Φ = Φ(ν)` (or of the lattice measure `μ`, with `ν` its KL scalar spectral measure) satisfying:

```
(Flip)  Φ(ν) ≤ inf supp ν = m₁²   for EVERY admissible (OS/RP) theory ν ,
        Φ ≡ 0 on the free massless field   (ν a.c. two-photon cut, f₁²=0, m₁=0),   [acid test]
        Φ ≥ λ > 0 uniformly in (L,a) on pure SU(N) YM ,
        and Φ is NOT (i) Schur-concave/moment-type, (ii) built from E₂/I2, (iii) deformation-invariant/topological/frozen.
```

The `E2 §6` version adds: "computed from trial data **richer than `(R,Var)`** without a certified near-eigenvector." The `HISTORICAL-TOOL-SWEEP §3` gives the *only* live shape it could not kill: **"a sector/threshold mismatch giving a genuinely kinematic `E₂`."** The present task replaces that with the **Savvidy/dynamical** shape. We attack that directly.

**Key structural discriminator the corpus already established.** The gap is **support** `m₁²=inf supp ν`, not the atom (`SUPPORT-NOT-ATOM`); equivalently the **spectral gap of the coarse-mode generator** `L = Δ − ∇V·∇` (Poincaré constant `ρ`), `m_eff²(k→0) ≥ ρ ⟺` H1. So a candidate `Φ` must lower-bound **the spectral gap of a diffusion generator**, not a Hessian bottom. The Savvidy object is the latter.

---

## §1 The candidate Φ — writing it down concretely

We make the most generous, most literal version of "use Savvidy as the lever." On the lattice measure `μ = Z⁻¹ e^{−V} dX` (`V = S_W + Φ_Haar`, `X = log U ∈ su(N)` on the gauge quotient `X = G^E/G_site`), let

```
   H_phys(x) := Π(x)† · Hess V(x) · Π(x)        (Hessian restricted to the gauge-invariant
                                                  transverse soft sector, `Π` the orbit quotient),
   Φ_Sav(μ)  := inf_{x ∈ supp μ}  λ_min( H_phys(x) )      = the deepest pointwise curvature
                                                             (the Savvidy depth, incl. Haar).
```

- On a constant abelian background `λ_min(H_phys) = −g|B|` (Landau-ladder bottom), the **Savvidy depth**.
- On the free field `F=0`: `Hess V = −Δ ⪰ 0`, `λ_min = 0 ⇒ Φ_Sav = 0` (acid-test letter ✔).
- On free massive Proca `–Δ+m²`: `λ_min = m²` (acid-test letter ✔).

**The attempted flip inequality (three natural forms — all that "lever" could mean):**

```
   (L1)   m₁²  ≥  −Φ_Sav            "the Savvidy depth is the mass scale"      (a lower bound via depth)
   (L2)   m₁²  ≥  c·|Φ_Sav| ,  c>0  "deeper Savvidy ⇒ bigger gap"              (monotone in the depth)
   (L3)   m₁²  ≥  max(0, Φ_Sav)     "if Hessian positive everywhere, gapped"    (positivity ⟹ gap)
```

(If none holds, one may still try `Φ :=` a *functional of the whole depth distribution* `dμ-width`, not just `inf`; §4 shows the anti-correlation defeats any monotone such map.)

**The decisive question is whether any (L·) is true on a non-convex (Savvidy-type) theory.** The corpus's own data already answer "no" for the real problem: the **static-curvature ("Savvidy") mass `~0.6β` is *not* the transfer-matrix gap** (`TEST-C-GLOBAL-SCHUR §4`; `TOOL-FACTORIZATION §4 Q3` explicitly excludes `~g|B|` as "not the transfer-matrix gap"). §2 proves this is not a coincidence of the real theory but a **structural theorem** any non-convex-action theory must obey.

---

## §2 The acid test — a solvable-model witness [T in-model] + [V-here]

To test `(L1–L3)` **exactly**, go to the simplest measure whose transfer gap is computable and whose action Hessian goes Savvidy-negative — a 1D non-convex "measure":

```
   μ = Z⁻¹ e^{−V(x)} dx ,   generator L = d²/dx² − V′(x) d/dx ,
   transfer gap  m₁²  := spectral gap of L on L²(μ)  = Poincaré constant ρ.
   Exact tool:  U(−L)U⁻¹ = H_s := −d² + ( V′²/4 − V″/2 )  on L²(dx)   [Uf = f e^{−V/2}],
                so  m₁² = E₁(H_s) − E₀(H_s) = E₁(H_s)   (E₀=0).
   "Savvidy depth"  Φ_Sav = inf_x V″(x) .
```

This is the **same** object the corpus uses for H1 (spectral gap of the coarse Langevin generator), reduced to a solvable instance; the "Hessian bottom" is the literal 1D Savvidy object. Code: `research/unconditional/e2_flip_toy.py` (headless, single process, `numpy`).

### (A) Free massless `V=0` — must give **no gap**
Hessian bottom `inf V″ = 0`. Dirichlet-box gap `= (2π/L_box)²−(π/L_box)² → 0` as `L_box → ∞` (computed `0.2056` at `L_box=12`; →0). **Both → 0. ✔ acid test passes (letter).**

### (B) Free massive `V=(1/2)m²x²`, `m=1` — must give **mass `m`**
Hessian bottom `= m² = 1.000000`. Measured transfer gap `= 1.000001 = m²` (`E₀=−3.2e−7`, `E₁=1.000000`). **`Φ_Hess = gap` EXACTLY. ✔ works perfectly on the Gaussian (`mass = 1`).**

### (C) Savvidy-like double well `V=(g/4)(x²−v²)²`, `v=1` — the decisive case
`V″(0) = −g v² < 0` (the Savvidy soft-mode, an exact negative Hessian at the symmetric point); the measure is **bimodal** (two "vacua" at `±v`); the transfer gap is the **inter-well tunnelling splitting** — the exact 1D image of the nonperturbative transmuation mass `~e^{−S_inst}`. The Hessian depth and the true gap (from `e2_flip_toy.py`, finer-grid confirmation in the deep regime):

| `g` | `Φ_Sav = inf V″ = −g` | **true gap `m₁²`** | mass `m₁` | `gap/|Φ_Sav|` |
|---|---|---|---|---|
| 0.25 | −0.25 | 5.313e−01 | 0.729 | 2.13e+00 |
| 1.0 | −1.0 | 7.921e−01 | 0.890 | 7.92e−01 |
| 4.0 | −4.0 | 7.487e−01 | 0.865 | 1.87e−01 |
| 8.0 | −8.0 | 4.545e−01 | 0.674 | 5.68e−02 |
| 16.0 | −16.0 | 1.191e−01 | 0.345 | 7.44e−03 |
| 24.0 | −24.0 | 2.485e−02 | 0.158 | 1.04e−03 |
| 36.0 | −36.0 | 1.906e−03 | 0.044 | 5.30e−05 |

**This kills `(L1)`, `(L2)`, `(L3)` simultaneously and decisively:**
- **`(L1)` `m₁² ≥ −Φ_Sav`:** on the free massive field it is an equality (`1 ≥ 1`); at `g=8` it reads `0.454 ≥ 8` — **FALSE by a factor 18**. The inequality is not "loose," it is **wrong-signed**.
- **`(L2)` `m₁² ≥ c|Φ_Sav|`:** `gap/|Φ_Sav| → 0` **monotonically** (`2.1 → 0.057 → 5e−5`); no `c>0` works. Deeper Savvidy gives **smaller** gap, the **opposite** monotonicity.
- **`(L3)` positivity ⟹ gap:** the reverse (negativity ⟹ no gap) is **false** — the gap is `>0` precisely where the Hessian is negative. And on the free massless field `(L3)` returns `0` while `Φ_Sav=0`, so it certifies nothing there (vacuity, the Family-I pattern).

**The exact identity check (Stein), same run:** `E_μ[H] = E_μ[V′²]` holds to 3–4 digits in every row (`g=1`: `2.125406` vs `2.125395`; `g=8`: `12.451607` vs `12.451197`). This **confirms the corpus's L3** (`H1-HESSIAN-DIRECT §1.1`: the unconditioned expected Hessian is PSD, `=E[∇V⊗∇V]`). **And `E_μ[H]=E_μ[V″]` is *positively large and increasing in `g`* (`≈2g`) while the true gap *decreases* — the expected Hessian is exactly as gap-blind in the solvable model as the corpus proved it is in QFT.** At the symmetric point the *conditioned* curvature is `G″(0)=V″(0)=−g<0` (the bimodal dip) while the true (global) gap is positive — the **L3/L4 conditioning separation** in a fully solvable model.

**Reading.** In the double well the mass is generated by **tunnelling** (`m₁ ≈ ω√(S_inst/2π) e^{−S_inst}`, `S_inst=(4v³/3)√(g/2) ∝ √g`), a **global/non-local** functional of the whole action landscape. The local second variation `V″` sees only the *barrier-top curvature* (negative) and the *well curvature* (positive); it is **anti-correlated** with the tunnelling rate. **This is the exact structural image of YM:** the Savvidy depth `−g|B|` grows with coupling while the physical glueball mass is exponentially small (`~e^{−c/g²}`) — a fact the corpus records numerically (`TEST-C`: `~0.6β` static ≠ transfer gap). **Savvidy negativity is the *cause* of the transmutation, hence necessarily *not* a lower bound on the resulting mass.**

---

## §3 Scorecard — the candidate `Φ_Sav` against the six flip constraints

| # | constraint (`E2 §6` + task) | `Φ_Sav` | basis |
|---|---|---|---|
| (i) | **not Schur-concave / not moment-type** | **PASS** | `Φ_Sav` is a genuine *nonlinear* functional of the measure (an inf of a curvature over `supp μ`); it is not in the linear-data cone. This is the one constraint it genuinely satisfies — *and it is why the task's suggestion is attractive.* |
| (ii) | **does not use `E₂`** | **PASS** | uses only `Hess V` (a zero-frequency datum), no second level. |
| (iii) | **coupling-dependent/dynamical (Savvidy)** | **PASS (but the wrong direction)** | `Φ_Sav=−g|B|` is β- and `g`-dependent; but §2 shows it moves **opposite** to the gap. |
| (iv) | **acid test (fails free massless)** | **PASS letter / VACUOUS** | free massless `Φ_Sav=0` (correct number) but it then certifies nothing (`0≤0`); free massive `Φ_Sav=m²=gap` — i.e. the *only* place it equals the gap is the trivial Gaussian. |
| (v) | **uniform in `(L,a)`** | **FAIL** | `N_neg ∝ V`, depth `∝ β`; `Φ_Sav` is extensive/UV (`PROOF-SAVVIDY-EXACT-SIGN`). §2's anti-correlation is already fatal before uniformity. |
| (vi) | **non-circular & direction-correct (`Φ ≤ m₁²`)** | **FAIL** | `(L1)` false by factor 18 at `g=8`; the direction-correct version (`Φ = e^{−S_inst}`) **is** the nonperturbative transport machine = H1. |

**Net `Φ_Sav`: 3/6 formal passes, 0 of them load-bearing.** It clears the *easy* three (i,ii,iii) — genuinely escaping the Schur/regress/topology walls — and dies on exactly the three that matter (iv-vacuity, v-uniformity, vi-direction). This is the sharpest possible confirmation of the `E2` diagnosis: **the constraints that the Savvidy lever escapes (Schur / E₂ / rigidity) are not the binding ones; the binding ones (acid-test non-vacuity + uniformity + direction) are precisely what the Hessian cannot supply.**

---

## §4 Why *every* Φ collapses — the exhaustion, and the precise incompatibility

The `E2 §6` flip condition and the task's six constraints, taken together, are inconsistent **for the Savvidy class**, and the corpus's four existing class-no-gos complete the exhaustion over the other named data-kinds:

1. **Moment/linear/transform functionals (Family I)** → Schur-concave, upper bounds only (`SCHUR-CONVEX`, `TOOL-COMPLETELY-MONOTONE`). Excluded by (vi)-direction. *The Savvidy lever's only advantage is escaping this — and it pays for it in (v)+(vi).*
2. **Variational/index/regress (Family II)** → need `E₂`/I2; in a one-scale theory `E₂=ρE₁` so "isolation" ≡ "gap" (`E2 §0.4`, `G2D-ISOLATION`). Excluded by (ii)/(vi).
3. **Deformation-invariant / topological / derived** → mass-blind, frozen at free value (`TOOL-FACTORIZATION`, `INVENT-NEW-INDEX` rigidity). Excluded by (iv)-non-vacuity + (iii).
4. **Raw action-Hessian / second-variation (the Savvidy lever, this doc)** → **anti-correlated with the gap** (witness §2), extensive (v), and its direction-correct form is the tunnelling action = H1 (vi). Excluded.

**The unifying structural statement (the real content of this document):**

> **[T, scoped]** *The E₂-flip functional, if it exists, must lower-bound the spectral gap of the **conditioned** measure — a global/non-local functional reflecting tunnelling between the unstable (Savvidy) configurations. The Savvidy Hessian is the **local second variation at a single configuration**; the two differ by exactly the **conditioning/variance term** `D = Var(f|0) − E[⟨K⁻¹c,c⟩|0]` of `G″ = S − D` (`INVENT-VARIANCE-INEQ`, `H1-HESSIAN-DIRECT §3`). The witness §2 is the solvable-model instance of this identity: `V″` (local) and `m₁²` (global) are separated by the tunnelling term, and control of that term **is H1**.*

So the Savvidy "lever" does not *fail to be tried well enough* — it **fails because the Hessian is, by construction, one derivative/one conditioning level away from the gap, and the missing level is exactly the open problem.** This is a genuine incompatibility of the constraint set (i,iii) vs (v,vi) for the named class, **not merely this attempt's weakness.**

**Is it `∀possible`? No — and we must be exact about that.** This is a **[NG] for the named class** "monotone/continuous functionals of `Hess V` (the action's local second variation) and their measure-averages," together with the four previously-banked class-no-gos. It is **not** a global barrier (`KEY-EXISTS §5`: `∀class ≠ ∀possible`). A *new* functional that lower-bounds the **conditioned** gap `G″=S−D` uniformly — i.e. that certifies `D<S−λ` — would be a genuine flip, but it is *by definition* H1 and no such functional is known (it is the exact stuck step, §5).

---

## §5 The exact stuck step / what would flip this

The single object whose construction flips the verdict (`=H1`, in three equivalent faces already in the corpus):

```
   (H1)   Uniform bound on the conditioned coarse force-variance:
            Var_μ( ∂_n V | a=0 ) − E_μ[H_vv | a=0]  <  W_full − λ      [= D < S − λ]
          ⟺  Var_μ(a_{k→0}) ≤ 1/λ                                    [Ward]
          ⟺  G″(k→0) = S − D ≥ λ > 0,   uniform (L,a), across β≈2.5.
```

A flip requires a functional of the **conditioned** measure — not of `Hess V` — that bounds `D = Var(f|0) − E[⟨K⁻¹c,c⟩|0]` from **above** without presupposing `Var(a)≤1/λ`. This is a **variance ceiling**, direction `≤`, on the *conditioned* slice. Every positive piece the corpus can prove (`W_geo≥0`, Haar floor `N/12`, the Stein identity) lives on the **unconditioned** side or on the **complement** of the sign-carrying cell (typical field × `k→0`), which is why `TRANSPORT-CORE §4`'s "anti-localization of positivity" bites. The Savvidy lever cannot reach it because the Savvidy datum is the *unconditioned* pointwise Hessian, on the wrong side of the conditioning that defines `D`.

**What would count as a genuine Savvidy-based flip (stated so it can be refuted):** a **provable inequality** `Var_μ(∂_n V|a=0) ≤ E_μ[H_vv|a=0] + W_full − λ` whose proof **uses the Savvidy spectrum** to bound the *conditioned* variance (e.g. by showing the Savvidy-negative directions are *decorrelated* in the conditioned slice so they inflate `S` more than `D`). The witness §2 says the local landscape alone cannot do this (there the negative direction *is* what makes `D` large / `G″(0)<0`). No such inequality is known; producing one is H1. **[O]**

---

## §6 Honest boundary — what this document does and does not say

- **Does NOT prove Clay; does NOT prove no tool can exist** (`∀class ≠ ∀possible`). It is a **[NG] for the named class** "functionals of the action's local Hessian / second variation" + the four banked class-no-gos, with the **fresh solvable-model witness [V-here]** that the Savvidy depth is *anti-correlated* with the gap (the reason "use Savvidy as a lever" is impossible, not merely unfruitful).
- **The witness is a solvable 1D model, labelled as such.** Its role is to test the *logic* of `(L1–L3)` exactly, not to compute the SU(N) gap. The structural content (`Hess V` ≠ transfer gap; gap = global tunnelling/conditioning; Stein `E[H]=E[V′²]` gap-blind) transfers dimension- and gauge-independently — it is the solvable instance of the corpus's already-established `TEST-C` / `H1-HESSIAN-DIRECT §4` facts. **No SU(N) number is claimed.**
- **The three constraints the candidate passes (i,ii,iii) are real and genuinely escape the Schur / E₂ / topology walls** — under-claiming them would be dishonest the other way. The verdict is that **the three it fails (acid-test non-vacuity, uniformity, direction) are the binding ones**, which is the `E2` thesis made concrete.
- **The "leak" — the direction-correct Savvidy functional `Φ = e^{−S_inst}` (the tunnelling/condensation action) — is a real object:** it is positive, dynamical, acid-test-clean, non-Schur, non-regress, and *equals* the gap. But it does not *flip* the problem: **computing it uniformly `(L,a)` across the Savvidy crossover is precisely the nonperturbative transport machine that does not yet exist** (`CLAUDE.md §5`, `TRANSPORT-MACHINE`, `W3`). So the Savvidy lever, done correctly, reduces to the missing machine — not past it.
- **No fabricated flip, no QED.** Clay stays open.

---

## References

**Internal (absolute paths; the verdicts this document classifies against and extends):**
- `D:\Dev\math\CLAUDE.md` (§1 build-don't-stop + acid-test iron rule; §2 labels; §5 wall & missing machine); `D:\Dev\math\docs\BRIEF-FOR-MODEL.md` (§2 five constraints; §3 death map; §7 "positivity proves disorder/upper; atom needs order/lower").
- `D:\Dev\math\research\unconditional\E2-REGRESS-CIRCULARITY.md` (§0.4 SSC; §0.6 two horns; **§6 the flip condition this doc attacks**); `KEY-EXISTS-OR-IMPOSSIBLE.md` (§5 `∀class≠∀possible`); `HISTORICAL-TOOL-SWEEP.md` (§3 the "kinematic-E₂" live shape; Family I/II); `TOOL-FACTORIZATION-ALGEBRA.md` (two-tier blindness; §4 Q3 `~g|B|` ≠ transfer gap); `TOOL-COMPLETELY-MONOTONE.md`.
- `LEVER-nonconvex-edge.md` (§0B linear-data lemma); `SCHUR-CONVEX.md` (master Schur no-go; free field `1≤0`); `G2D-ISOLATION-I2.md` (SSC `E₂=ρE₁`).
- `PROOF-SAVVIDY-EXACT-SIGN.md` (`H=−D²+2adF`, ladder `−g|B|`, `N_neg∝V`, `Tr H≥0`); `H1-HESSIAN-DIRECT.md` (**Stein `E[H]=E[f²]⪰0` gap-blind; gap = L3→L4 conditioning `G″=E[H]−Var+W`**; `D/S≈2%`); `TRANSPORT-CORE.md` (`Γ₂≡H1` lossless; §4 anti-localization of positivity); `PROOF-JACOBIAN-SAVES.md` (`W_geo≥0`, `G″=Δ_dyn+W+X`); `INVENT-VARIANCE-INEQ.md` (`G″=S−D`); `TEST-C-GLOBAL-SCHUR.md` (§4 static-curvature `~0.6β` ≠ transfer gap); `PROOF-P3-UNCONDITIONAL.md` (strong-coupling residue `f₁²>0`, `β<1/48` — the only regime where the gap is a theorem).
- `new-tools-registry.md` (A2 `G″=S−D`; A5 `Γ₂`/W-entropy; B5 Krein; D-residual dictionary: H1, (CP), `f₁²`, `G″`, Savvidy sign); `e2_flip_toy.py` (**the [V-here] solvable-model witness**).

**External — [R] standard / textbook (not re-verified this session):**
- Savvidy, *Phys. Lett. B* **71** (1977) 133; Nielsen–Olesen, *Nucl. Phys. B* **144** (1978) 376 (chromomagnetic instability, Landau-ladder `−g|B|`).
- Bakry–Émery, *Diffusions hypercontractives*, LNM 1123 (1985); Bakry–Gentil–Ledoux, *Analysis and Geometry of Markov Diffusion Operators*, Grundlehren **348** (2014) (Poincaré gap of a diffusion generator = spectral gap of the ground-state-transformed Schrödinger operator; the `U(−L)U⁻¹` transform used in §2).
- Coleman, *Aspects of Symmetry* (1985) / Callan–Coleman (1977); Garg, *Am. J. Phys.* **68** (2000) 430 — instanton/tunnelling splitting `ΔE ~ e^{−S_inst}` for the symmetric double well (the §2(C) gap).
- Prékopa, *Acta Sci. Math.* **32/34** (1971/73); Brascamp–Lieb, *J. Funct. Anal.* **22** (1976) 366; Stein's lemma / Schwinger–Dyson IBP (the `E[H]=E[∇V⊗∇V]` identity).
- Shen–Zhu–Zhu, *CMP* **400** (2023) 805, arXiv:2204.12737 (4D lattice YM gap `β<1/48`).

---

*Status: **[NG] for the named class** "monotone/continuous functionals of the action's local second variation (`Hess V`) and their measure-averages" (Savvidy lever), on the strength of a **[T]-in-model + [V-here]** anti-correlation witness (`Hess` depth `→−∞` while the true transfer gap `→0` exponentially), plus the four banked class-no-gos; **[O] exact missing step** = the conditioned force-variance ceiling `D<S−λ` `=H1` (uniform, across `β≈2.5`); **[J]** that (i,iii) vs (v,vi) is a genuine incompatibility for this class. Verdict: **IRREDUCIBLE — no flip constructed; the Savvidy "lever" is provably the wrong direction, and the direction-correct Savvidy functional (`e^{−S_inst}`) is the missing transport machine under its standard name.** NOT a proof of Clay; NOT a global impossibility theorem. **Clay stays open.***
