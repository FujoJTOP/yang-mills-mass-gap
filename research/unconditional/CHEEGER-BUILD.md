# CHEEGER-BUILD — Cheeger/Buser isoperimetric spectral gap for 4D pure SU(N) YM: the one species whose LOWER bound never signs the indefinite Hessian (real structural gain), but whose constant `h>0` uniform is *verbatim* H1, because Cheeger's infimum is a **global** object attained by the k→0 softest-mode half-space — local tube isoperimetry survives the Savvidy cone (the cone is mass-starved, WSTC), the global minimizer does not.

> **2026-11-23 · CONSTRUCTOR build. 判词 (verdict).** Cheeger is the only species in the family whose gap bound `λ₁ ≥ h²/4` is **sign-free** — it never asks `G″⪰0`, never signs `−D²+2adF_B`, and therefore does *not* die at the Savvidy wall the way Bakry–Émery (`Γ₂<0`) and Steklov (indefinite DtN) do. That advantage is real and I certify it: on the classical Wilson measure the Savvidy negative mode lives on the mass-starved cone (`ν(C_κ)≤e^{−cβΔκ²N_□}`, WSTC), and a **thin-tube concentration–Cheeger (CTC)** inequality passes all its hypotheses at β=2.5, returning positive *local* isoperimetry. **But** Cheeger's constant is a **global infimum** `h=inf_A μ⁺(∂A)/min(μ(A),μ(Aᶜ))`; its minimizer is not a tube or the cone — it is the k→0 **softest-mode half-space** `{ã(0)<c}`, whose ratio is `≍1/√(Var_μ ã(0)) = 1/√χ`. Hence `h(μ)≥h₀>0` uniformly in `(L,a)` **⟺** `χ≤1/h₀²` uniformly **⟺** H1. Cheeger dodges the curvature wall and lands on the bottleneck wall, and the bottleneck *is* the gap. **Passes #1, #4; #2 partial (local sign-free lower bound, does not become global); fails #3, #5. Score 2/5.**

---

## §0 · Step-0: what is already banked, and what is genuinely new here

`Grep` for `Cheeger / isoperimetr / Buser / Poincaré / Dirichlet form / spectral gap` returns 100+ files. The **Cheeger species is not virgin** — five banked results already fence it, and this build *reuses* them, does not re-derive:

| banked file | its result | relation to this build |
|---|---|---|
| **B3-ISOPERIMETRIC** | Cheeger on the *marginal* `ν`: fiber-saturated competitors carry **zero** fiber boundary ⇒ `h(ν)` = transverse bottleneck = H1. 2/5. | **the marginal half** of the story; I complete it with the *full-link* metric and the global-inf argument. |
| **POTENTIAL-THEORY-CAPACITY** | equilibrium/Wiener capacity on the 0-mode marginal; finite-volume vacuity; Arrhenius barrier `H≈0.3/λ=O(1)`; **non-vacuous** free-field acid test. 2/5. | the **free-field discriminator** and the Arrhenius shell; adopted verbatim. |
| **C-BBS-FORMULATE / D-BBS-MEASURE-TRAP** | capacity/conductance **trap**: cone killed by mass (WSTC), thermal shell is the only bottleneck (`≈0.55` per-cut at one scale), `k→0` uniformization = H1. | the **two-region conductance** picture; adopted. |
| **LOC-ABSOLUTE-CRUX** | fixed-cell Cheeger/Sinclair is **acid-test-blind** without L-uniformity (`V̄≡0⇒H₀=0` passes it); with L-uniformity it *is* H1. No seedless reading. | the **global-vs-fixed-cell** no-go; adopted. |
| **PROOF-CONCENTRATION** | concentration is **downstream** of the gap (Herbst one-way); volume, not curvature, kills uniformity. | the **concentration-escape** no-go (§2 states why CTC survives *locally* but not *globally*). |

**Genuinely new in this build:** (i) the **explicit Cheeger/Buser statement on the full link space** with the link-metric perimeter and the isoperimetric profile (§1); (ii) the derivation that the **global Cheeger infimum is governed by the k→0 mode** with the precise ratio `h(μ)≍1/√χ` (§1, the load-bearing computation); (iii) the **thin-tube concentration–Cheeger (CTC)** candidate as a precise conditional inequality, and the demonstration that its **local** positivity is genuine but does not reach the **global** minimizer (§2); (iv) the three-regime table and the scorecard (§3–§4). Everything else is a faithful synthesis of the five banked verdicts, cross-referenced.

---

## §1 · Cheeger/Buser for the YM generator, explicitly

### §1.1 Space, metric, Dirichlet form, generator
Configuration space `X = G^{E}`, `G=SU(N)`, `E = d·L^d = 4V` oriented links on `T^4` (`V=L^4` sites; each link a group element `U_ℓ∈G`). `dim X = E·(N²−1)`.

**Link metric.** Product of bi-invariant metrics, `d(U,U')² = ∑_{ℓ∈E} d_G(U_ℓ,U'_ℓ)²`, `d_G(e^{iθ})=|θ|` on the maximal torus extended bi-invariantly (`−2 tr(U_ℓ^†U'_ℓ)/N` up to normalization). `X` is compact ⇒ complete, finite volume; Cheeger's inequality applies with **no completeness caveat**.

**Measure.** `dμ = Z^{−1}e^{−S}dU`, `dU` = product Haar on `G^E`, `S = (β/4)∑_{□} tr(1−Re U_□)` on the lattice (`≅(1/2g²)∫trF²` in the continuum; `β=2N/g²`).

**Dirichlet form / generator.** Let `D_ℓ` be an orthonormal frame of right-invariant derivatives on link `ℓ`. The reversible (Glauber/heat-bath) Dirichlet form and generator are
```
ℰ_μ(f) = ∑_{ℓ∈E} ∫_X |D_ℓ f|² dμ ,        L = −∑_{ℓ} D_ℓ^* D_ℓ ,     μ-reversible,
λ₁(μ)  = inf_{Var_μ f>0} ℰ_μ(f)/Var_μ(f) .
```
`λ₁(μ)` is the **full-link dynamical gap** (slowest relaxation of the link heat-bath). It is the object Cheeger lower-bounds. (H1's `λ` is the *coarse-marginal* gap `λ₁(ν)`; §1.4 shows they agree on the softest mode, which is why the distinction is not an escape.)

### §1.2 Cheeger (1970), discrete form (Jerrum–Sinclair 1989, Sinclair 1992)
```
λ₁(μ) ≥ h(μ)²/4 ,      h(μ) = inf_{A⊂X} μ⁺(∂A) / min(μ(A), μ(Aᶜ)) ,
μ⁺(∂A) = ∑_{ℓ} ∫_X |D_ℓ 1_A| dμ        (link-metric perimeter; coarea),
I_μ(v) = inf{ μ⁺(∂A) : μ(A)=v }   (isoperimetric profile),   h(μ)=inf_{v∈(0,1)} I_μ(v)/min(v,1−v).
```
**What `μ(∂A)` means here.** `∂A` is a codim-1 hypersurface in `X`; its measure `μ⁺(∂A)` is the **total variation of `1_A` under the link-metric gradient** — the expected (over `μ`) total size of the jump in `1_A` when each link is moved a unit step, summed over all links. Equivalently, by coarea, `μ⁺(∂A)=∫_0^1 (perimeter density of the level sets of 1_A) dv`. It is the boundary of `A` in the *link* geometry, not in the coarse geometry. **Crucially: no curvature, no Hessian, no sign of `−D²+2adF_B` enters Cheeger's RHS** — `|D_ℓ 1_A|²⩾0` for any action.

### §1.3 Buser (1982) — the converse needs Ricci, i.e. the very object Savvidy breaks
```
λ₁(μ) ≤ 2√(δ(n−1))·h(μ) + 10 h(μ)²     provided   Ric_g ≥ −(n−1)δ ,  n=dim X.
```
The measure-reweighted Ricci of `(X, link metric)` w.r.t. `μ` is `Ric + Hess S`. At β=2.5 `min(Ric+Hess S_W) ≈ −2.3` (W6-BE; `Hess` reaches `−1.12β`), so the Buser hypothesis holds only with `δ≈2.3`, and the **converse direction is dead at the crossover**. This is *in our favour*: the direction we need — `h>0 ⇒ λ₁>0` — is **Cheeger**, which is unconditional. The asymmetry `Cheeger unconditional / Buser conditional-on-Ric` is the structural reason the species is worth building at all: **it survives `Ric<0` where BE dies.** The catch is below (§1.4): what it survives *into* is H1.

### §1.4 The load-bearing computation — the global inf is the k→0 mode, so `h(μ)≍1/√χ`
The k→0 coarse observable is the Fourier zero mode `ã(0) = V^{−1/2}∑_x a(x)` (`a` a component of the link/site field), with `Var_μ(ã(0)) = χ` (the susceptibility; `χ=1/m²` for a free massive field, `χ=∞` massless). H1 is `χ ≤ 1/λ`, uniform `(L,a)`.

- **Dirichlet energy of the mode is `O(1)`, not diluted.** Changing one link affects `O(1)` sites of `∑_x a(x)`, so `|D_ℓ ã(0)| = O(V^{−1/2})` per link; summing over `E=4V` links, `|∇ã(0)|² = ∑_ℓ |D_ℓ ã(0)|² = O(1)`. (Coarse summability is exact: the gradient is spread as `V` links × `V^{−1}`, total `O(1)`.)
- **Its level sets are the global Cheeger competitors.** For `A_c = {ã(0)<c}`, coarea gives
  `μ⁺(∂A_c) = ∫ δ(ã(0)−c)|∇ã(0)| dμ = p(c)·E[|∇ã(0)| | ã(0)=c] = p(c)·O(1)`, with `p` the marginal density of `ã(0)`. At the median `min(μ(A_c),μ(A_cᶜ))=½`, so
```
        μ⁺(∂A_c)/min(μ(A_c),μ(A_cᶜ))  =  2 p(median)·O(1)  ≍  1/σ_ã(0)  =  1/√χ .
```
- Since `h` is the **infimum** over *all* `A`, `h(μ) ≤ inf_c (ratio) ≍ 1/√χ`. Conversely, for a log-concave (or near-Gaussian) marginal the half-space is the minimizer, so `h(μ) ≍ 1/√χ` up to absolute constants.
```
        h(μ) ≥ h₀ > 0  uniformly (L,a)   ⟺   χ ≤ C²/h₀²  uniformly   ⟺   H1 .
```
**This is the entire content of the species.** Cheeger's sign-blindness is real (§1.3), but the object whose positivity it needs, `h(μ)>0`, is governed by the **softest mode's variance `χ`**, i.e. it *is* H1 in geometric costume. The full-link gap and the coarse gap coincide on this mode; there is no metric escape. This sharpens B3-ISOPERIMETRIC (marginal only) to the full space, and it is the reason every later section terminates at H1.

---

## §2 · The central attempt — the concentration escape, made precise (CTC)

**Goal of the attempt.** Cheeger survives `Ric<0`; the Savvidy cone is exactly `Ric<0` with the measure concentrated *away* from the bad directions (`S_W=‖F‖²` is minimized at `F=0`, `Hess S_W=−D²⪰0` there; the negative `−g|B|` is a one-loop/effective/dynamical artifact at nonzero background, i.e. on the cone where `μ` is starved). So: **can the concentration of `e^{−β‖F‖²V}` on near-minimal `|F|` deliver `h>0` despite the indefinite Hessian?**

### §2.1 The candidate inequality (thin-tube concentration–Cheeger, CTC)
Let `M = {F=0}` ⊂ `X` be the flat-connection set (a submanifold; on `T^4`, dimension `4(N−1)` modulo gauge), `M^r = {dist(A,M)≤r}`. Hypotheses:
- **(H-a) tubular concentration:** `μ(X∖M^r) ≤ ε_r`, with `ε_r` super-polynomially small in `β`.
- **(H-b) transverse positivity on the measure-carrying region:** `Hess S_W ⪰ κ·g` on the normal bundle of `M^r`, `κ=κ(β)>0`.
- **(H-c) intrinsic positivity:** the induced (conditional) measure `μ(·|M^r)` has Cheeger constant `≥ h_M>0` w.r.t. the restricted link metric.

**CTC (candidate, conditional):**
```
     h(μ)  ≥  ½·min(h_M, √κ)·(1 − 2ε_r)  −  C·ε_r^{1/2}/r ,
     hence  λ₁(μ) ≥ h²/4 > 0   provided  ε_r<½,  h_M,κ bounded below uniformly in (L,a).
```
(The constants are schematic — the content is the hypothesis-check and the identification of *which* part of `h` this controls. The tube/first-variation heuristic: a measure on a tube with normal profile of width `1/√κ` and escaped mass `ε_r` has Cheeger constant governed by the smaller of the intrinsic `h_M` and the normal width `√κ`, up to the `ε_r` leak.)

### §2.2 Hypothesis-check at β=2.5 — and what CTC does and does not reach
- **(H-a) HOLDS (banked).** WSTC (ABC-NEW-INEQ) proves the high-field cone `ν(C_κ)≤e^{−cβ(κ²−κ_typ²)N_□}`, uniformly in `L`. The large-field tail is super-polynomially starved, so `ε_r` is genuine: **`μ` concentrates near `M`.** This is the concentration the task asks for, and it is a theorem, not a heuristic.
- **(H-b) HOLDS on the measure-carrying region.** On `M^r` (`F≈0`) the classical Hessian is `−D²+2adF_B⪰−D²⪰0` (transverse modes gapped at `~1/a`); the Savvidy negative `−g|B|` sits on the cone (`r` large), where `μ` has mass `≤e^{−cβN_□}` (WSTC). So the indefinite direction is **at the mass-starved rim, not in the bulk** — precisely why Cheeger (sign-blind) is unaffected where BE (which weights by `e^{−S}` and signs `Hess S`) is killed. **This is the real gain: "isoperimetry can be positive under negative curvature" is TRUE locally, and CTC certifies it.**
- **(H-c) HOLDS formally.** The flat-connection moduli carry a compact flat measure; `h_M>0` (a compact flat torus has positive Cheeger constant `~1/diam`).

**So CTC passes all three hypotheses at β=2.5 and yields positive *local* isoperimetry.** It genuinely neutralizes the Savvidy cone: the cone is the wrong enemy, killed by mass (`WSTC`), not by curvature; the surplus convexity of the Jacobian (`W_geo=E[Φ″]−Var Φ′≥0`, Prékopa, PROOF-JACOBIAN-SAVES; pointwise floor `Φ″⪰N/12`) reinforces the bulk. **This part is a conditional theorem** (given H-a,H-b,H-c, which hold).

### §2.3 Why CTC does not close — the global infimum is not in the tube
Cheeger's `h` is `inf_A` over **all** measurable `A`, and the infimum is _not_ attained by local tube sets or the cone. By §1.4 it is attained by the **global k→0 half-space** `A_c={ã(0)<c}`, whose perimeter/min-mass ratio is `≍1/√χ`. That half-space is a **long-wavelength hypersurface spanning the whole torus** — it is *not* localized in `M^r`, so CTC's local `h_M` and `√κ` say nothing about it. Formally: CTC bounds the Cheeger ratio on the subclass of sets supported in `M^r`; the global inf can be (and is) smaller:
```
     h(μ) = inf_{A⊂X} ratio  ≤  inf_{c} ratio({ã(0)<c})  ≍  1/√χ .
     CTC controls  inf_{A⊂M^r} ratio  ≥ min(h_M,√κ) ,   but  min(h_M,√κ)  ≱  1/√χ  in general,
     and at the target 1/√χ is the smaller (softest-mode) object.
```
**The escape fails at the global-vs-local gap in the Cheeger infimum, and the global minimizer is the marginal/thermal-shell bottleneck = H1.** Concretely at β=2.5: the thermal shell `|a|~λ^{−1/2}` is a **bulk saddle** `V̄″≈−0.3`; its per-cut conductance is the finite Arrhenius factor `e^{−2H}≈0.55` at *one* scale (C-BBS-FORMULATE, D-BBS-MEASURE-TRAP) — `O(1)`, not exponential, and the cone does **not** contribute to this cut (WSTC). But `0.55` is a **one-scale** conductance; making it an `(L,a)`-uniform global floor is exactly `Var(a)≤1/h₀²`, i.e. H1 (LOC-ABSOLUTE-CRUX: the fixed-cell floor is acid-test-blind without L-uniformity, and `V̄≡0` passes it). **The concentration is real, the cone is neutralized, and H1 is relocated — not removed — to the global (k→0) uniformization of the thermal-shell conductance.**

---

## §3 · Three-regime test (honest, exact numbers)

| regime | measure / structure | Cheeger data | verdict |
|---|---|---|---|
| **STRONG COUPLING `β<1/48`** | SZZ LSI(λ_μ) uniform in `V` (CMP 400 (2023) 805, arXiv:2204.12737); marginal `ν` near log-concave; confining scale `O(1)` in lattice units | `h(ν) ≳ √λ_μ > 0`, `λ₁ ≥ h²/4`; matches **Lemma SC** `λ=C_blk^{−2}λ_μ` and **P3** `f₁²>0` (PROOF-P3-UNCONDITIONAL) in sign and scaling | **AGREES with P3/SC.** But *consumes* SZZ's LSI as input — a faithful restatement, **not an independent proof**. (Same functional `f₁²`, same window.) |
| **FREE MASSLESS `g=0`, `Φ=0`** | `μ` = massless Gaussian; `V̄≡0`; k→0 mode flat | `h(μ)=0` (`σ=∞`), `Cap=0`, `Var=∞`, `λ₁=0` | **`h=0`, gap 0 — correctly NO positivity.** Acid test **met**, and in the capacity/CTC form it is **non-vacuous** (the flat mode has literally zero capacity; cf. POTENTIAL-THEORY-CAPACITY §4 — the cleanest discriminator in the species family). |
| **FREE MASSIVE (mass `m`)** | massive Gaussian; marginal `N(0,1/m²)` | `h ≍ √(2/π)·m`, `λ₁=m²` (marginal normalization); correlator decays `e^{−m|x|}` | **`h>0`, gap `=m²` (= mass `m` in the correlator).** Cheeger passes the massive control cleanly and uniformly in `L`. |
| **CROSSOVER `β≈2.5` (hard case)** | `PᵀHP` min `≈−0.3` (MARGINAL); `G″=+0.51` SU3 (W7 Part B); `D/S≈2%` (SCBL); `Ric+HessS≈−2.3`; cone mass `≤e^{−cβΔκ²N_□}` (WSTC); shell saddle `V̄″≈−0.3`, `e^{−2H}≈0.55` (one scale); `c<1` false, `ρ≈3–16` (PROOF-SAVVIDY-WGEO-DOMINATION) | **Local (CTC):** `h_local ≥ min(h_M,√κ)>0` — **positive**, cone neutralized by WSTC mass, Savvidy direction mass-starved. **Global:** `h(μ) ≤ 1/√χ`, minimizer = k→0 thermal-shell half-space; `0.55` is one-scale, **not L-uniform** | **LOCAL isoperimetry SURVIVES (real gain over BE); GLOBAL `h` collapses to `1/√χ` = H1.** The negative-mode direction does **not** collapse the *local* boundary fraction (it is on the mass-starved cone), but it does drain the *global* conductance into the thermal shell, whose uniform floor is H1. |

**Exact crossover statement.** The Savvidy negative mode is at the **rim** of the measure (`dist to cone large ⇒ mass `≤e^{−cβN_□}}`), so it does **not** set `h`; the **bulk thermal shell** does, and the shell's conductance is a bulk saddle, `H=O(1)` at `λ=O(1)`. So the negative mode `V̄″=−2adF_B` does **not** collapse the isoperimetric constant directly — the isoperimetric constant collapses instead through the **softest-mode variance** `1/√χ`, which is H1. **This is the precise difference between Cheeger and BE: BE is killed at the cone by the Hessian sign; Cheeger is killed at the k→0 marginal by the variance, with the cone removed from the problem entirely (WSTC). Wrong enemy for BE, right enemy (the marginal) for Cheeger — and the marginal is H1.**

---

## §4 · Five-condition scorecard (honest)

| # | condition | verdict | reason |
|---|---|---|---|
| 1 | **Nonperturbative** | **PASS** | Cheeger/capacity is exact, all `β`, no loop expansion; `|D_ℓ 1_A|²⩾0` for any action. |
| 2 | **Sign-controlling LOWER bound** | **PARTIAL** | Genuinely sign-free — the *only* species that lower-bounds the gap without signing `−D²+2adF_B`; survives `Ric<0` and certifies **local** positivity past the Savvidy wall (CTC §2.2). But on compact finite volume the bound is trivial (`h>0` vacuous), and the *useful* `h₀>0` is the global marginal = #5. Real operator-level pass, value-level fail. |
| 3 | **Non-circular** | **FAIL** | The uniform `h₀>0` **is** H1 (§1.4); the only proven regime (β<1/48) *consumes* SZZ's LSI (the gap); concentration is downstream of the gap (PROOF-CONCENTRATION). A faithful rename of H1, zero slack. |
| 4 | **Non-SUSY + fails free massless** | **PASS (cleanest)** | Pure YM, no superpartners; free massless `V̄≡0 ⇒ h=0, Cap=0, Var=∞` — fires on the free field, **non-vacuously** (capacity form). Discriminates YM from free. |
| 5 | **(L,a)-uniform `m_eff²(k→0)≥λ`** | **FAIL** | `h₀>0` uniform ⟺ `χ≤1/h₀²` uniform ⟺ H1; the thermal-shell conductance is `O(0.55)` at one scale, collapses in the `k→0`/continuum limit. |

**Score: 2/5** (strict; #1, #4 pass; #2 is a real but local sign-free lower bound that does not become the global uniform bound, so it rounds to a fail for closing). This is consistent with the banked B3-ISOPERIMETRIC (2/5) and POTENTIAL-THEORY-CAPACITY (2/5), and *improves on* them by isolating the exact structural asset (#2 local sign-blindness) and the exact cost (the global-vs-local gap in the Cheeger infimum).

---

## §5 · Sharpest line, exact stuck step, refs

**Sharpest line.** *Cheeger is the one species that never signs the indefinite Hessian, and it makes good on that: the Savvidy cone is killed by mass (`WSTC`), the bulk is kept positive by the Jacobian (`W_geo≥0`), and a thin-tube concentration inequality passes every hypothesis at β=2.5 and returns genuine local isoperimetry. But `h` is a global infimum, not a local one, and its minimizer is the k→0 softest-mode half-space — so `h>0` uniform is `χ≤1/h₀²` uniform is H1. Cheeger trades "prove `G″⪰0`" (curvature, false on the cone) for "prove `h>0`" (no bottleneck), and the no-bottleneck statement is the same gap. Wrong wall dodged, right wall hit: the bottleneck is the mass.*

**Exact stuck step.**
```
     h(μ) ≥ h₀ > 0  uniformly in (L,a)
  ⟺  inf_A μ⁺(∂A)/min(μ(A),μ(Aᶜ)) ≥ h₀  uniformly
  ⟺  the global minimizer (k→0 thermal-shell half-space {ã(0)<c}) has ratio ≥ h₀
  ⟺  1/√χ ≥ h₀  ⟺  Var_μ(ã(0)) ≤ 1/λ,  λ=h₀²/4   ⟺   H1.
Local (CTC) half proven/conditional:  inf_{A⊂M^r} ratio ≥ min(h_M,√κ) > 0  [cone WSTC-starved, bulk W_geo≥0].
Missing:  the global infimum  inf_{A⊂X} ratio  equals the marginal/thermal-shell value 1/√χ  —  i.e. H1 again.
```

**Single regime/step where it first breaks.** The **crossover β≈2.5, at the global infimum** (equivalently C-BBS/D-BBS `(R2)` no-thin-neck + `(R4)` k→0 uniformity): the *local* tube isoperimetry is positive and survives the Savvidy cone, but the *global* Cheeger constant is `1/√χ`, the softest-mode variance, which is exactly H1. (Strong coupling agrees with P3/SC; both free-field controls pass cleanly.)

**References.**
- *External (verified):* J. Cheeger, "A lower bound for the smallest eigenvalue of the Laplacian," in *Problems in Analysis* (Symp. in Honor of S. Bochner), Princeton Univ. Press (1970) 195–199. P. Buser, "A note on the isoperimetric constant," *Ann. Sci. École Norm. Sup.* **15** (1982) 213–230. A. Sinclair, "Improved bounds for mixing rates of Markov chains and multicommodity flow," *Combin. Probab. Comput.* **1** (1992) 351–370. M. Jerrum, A. Sinclair, *SIAM J. Comput.* **18** (1989) 1149. P. G. Doyle, J. L. Snell, *Random Walks and Electric Networks*, MAA Carus Math. Monogr. **22** (1984). M. Ledoux, *The Concentration of Measure Phenomenon*, AMS Math. Surveys Monogr. **89** (2001). H. Eyring, *J. Chem. Phys.* **3** (1935) 107; H. A. Kramers, *Physica* **7** (1940) 284. A. Prékopa, *Acta Sci. Math. Szeged* **32** (1971) 301; H. J. Brascamp, E. H. Lieb, *J. Funct. Anal.* **22** (1976) 366. R. Holley, D. Stroock, *Ann. Probab.* **15** (1987) 1541. **Shen–Zhu–Zhu**, *Comm. Math. Phys.* **400** (2023) 805, arXiv:2204.12737 (β<1/48 LSI). Savvidy, *Phys. Lett. B* **71** (1977) 133.
- *Internal (banked, cross-referenced):* **B3-ISOPERIMETRIC** (marginal Cheeger = H1, fiber-saturated competitors zero fiber boundary); **POTENTIAL-THEORY-CAPACITY** (equilibrium/Wiener capacity, finite-volume vacuity, non-vacuous free-field test, Arrhenius `0.3/λ`, `λ₁≥Φ_eq²/4`); **C-BBS-FORMULATE** (capacity trap, `E_ν` sign-free, cone mass-starved, thermal-shell `≈0.55`); **D-BBS-MEASURE-TRAP** (two-region trap invariant; gap-invariance = k→0 uniformization = H1); **LOC-ABSOLUTE-CRUX** (fixed-cell Cheeger acid-test-blind without L-uniformity, `=H1` with it); **PROOF-CONCENTRATION** (Herbst one-way; volume kills uniformity, not curvature); **ABC-NEW-INEQ** (WSTC cone starvation `e^{−cβΔκ²N_□}`); **PROOF-JACOBIAN-SAVES** (`W_geo=E[Φ″]−Var Φ′≥0`, `Φ″⪰N/12`, `G″=Δ_dyn+W_geo+X`); **PROOF-SAVVIDY-WGEO-DOMINATION** (`c<1` false, `ρ≈3–16`, log-cut anti-correlation); **OTHER-H1-POINCARE** (coarse form = H1, `PᵀHP` min −0.3); **PROOF-COARSE-POINCARE** (Lemma SC β<1/48; `G″` exact); **MARGINAL** (`PᵀHP` min −0.3); **W7-Gpp0** (`G″=+0.51`, SU3 β=2.5); **W6-BE-at-2.5** (`Ric+Hess S≈−2.3`); **PROOF-P3-UNCONDITIONAL** (`f₁²>0` at β<1/48); **INVENT-VARIANCE-INEQ** (SCBL `G″=S−D`, `D/S≈2%`).
