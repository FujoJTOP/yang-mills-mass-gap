# TRANSPORT-BLIND — independent blind red-team on the "Bakry–Émery Γ₂-calculus / optimal-transport (Otto–Wasserstein)" species for the 4D pure SU(N) mass gap

> 2026-11-23 · **判词：SPECIES DEAD. All four kill targets land; survival score 1/5.** The proposed tool claims the Γ₂ criterion (`Γ₂≥ρ ⇒ spectral gap ≥ρ`) is "the one positivity machinery that natively gives a gap **lower** bound," and that `Γ₂ = Ric_Haar + Hess_S` splits into a **proven-nonnegative Haar part** (`W_geo≥0`, Prékopa) plus a **Savvidy-indefinite action part** that "may be rescued by measure-starvation (WSTC)." It dies for the reason every positivity tool before it died: on the flat coarse space where the gap lives, `Γ₂` is **lossless** and equals `G″` **exactly** — so "`Γ₂≥ρ`" *is* "`G″⪰ρ`" *is* `m_eff²(k→0)≥ρ` *is* **H1**, with zero slack (B3-FUNCTIONAL-BE §1). The proposed decomposition is not a new lever — it is the exact `G″=Δ_dyn+W_geo+X` identity already banked (PROOF-JACOBIAN-SAVES), dressed in Riemannian geometry. The "rescue" (WSTC) provably starves the **wrong tail**: it suppresses the high-field cone, where the negative weight is measure-suppressed anyway, and leaves untouched the sign-deciding negative variance at **typical field, high momentum** (ABC-NEW-INEQ §4, W6-BE §3). And the "proven-nonnegative Haar part" is precisely the compactness artifact (`Ric` of SU(N)) that is a **kinematic** mass, not the dynamical gap, and vanishes in the continuum (B3 §5 caveat, W6 §1). This species is a faithful rename of H1, not a lowering of it. **Kills: (1) BE≡Hessian wall ✓; (2) starvation-on-wrong-tail ✓; (3) free-field/vacuity leak ✓; (4) circularity ✓.** Exact stuck step: `Var_μ(⟨v,∇S⟩|a=0) ≤ E[V″|a=0]−ρ` at typical field — (◆) of PROOF-JACOBIAN-SAVES, = H1.

---

## §0. The tool, stated precisely (kill the real thing, not a strawman)

Two faces of one species, both already probed in-repo and both false — the proposal is their union:

- **Full-space BE** (W6-BE-at-2.5): on the gauge quotient `M=G^E/G_site` with bi-invariant metric, `μ_Λ=Z^{−1}e^{−S_W}dU`, the claim is `CD(m,∞)` uniformly in volume, i.e. `Ric_M(v,v)+Hess S_W(v,v) ≥ m|v|²` ⇒ gap `≥ m` (Otto 2001; Sturm 2006; Lott–Villani 2009; von Renesse–Sturm 2005) ⇒ volume-uniform gap via LSI (Shen–Zhu–Zhu, CMP 400 (2023) 805).
- **Coarse/conditioned BE** (B3-FUNCTIONAL-BE): on `ν=P_*μ` (the block-averaged pushforward), the exact coarse generator `L_ν=Δ_c−∇S_eff·∇_c`; the criterion `Γ₂≥ρ` for the coarse mode `a_{k→0}`.
- **The "Perelman W-entropy" face** (INVENT-NEW-MONOTONE): the same `Γ₂=Ric_A+Hess_A S_W` as the (only) gap-certifying monotone.

The proposal's *novel* wrinkle is the **split + rescue**: `Γ₂ = Ric_Haar + Hess_S`, with `Ric_Haar` (= `W_geo≥0`, Prékopa) "proven nonnegative" and `Hess_S` (= Savvidy, indefinite) "rescued by WSTC measure-starvation." To be a *tool* it must (1) be nonperturbative; (2) deliver a sign-controlling **lower** bound `m_eff²(k→0)≥ρ>0`; (3) be non-circular (not assume the gap/`m_g>0`); (4) fail on the free field (`f₁²=0`, massless); (5) be `(L,a)`-uniform and reach `β≈2.5`. §1–§4 are the four kills; each maps to a hard constraint.

---

## §1. KILL 1 — Bakry–Émery is the Hessian wall renamed (`Γ₂=Ric+Hess_S` **is** `Hess⪰ρ`)

**Lands — decisive, on both faces.**

**(a) Full space = the Savvidy-sign wall, verbatim.** W6-BE-at-2.5 tests exactly `Ric_M(v,v)+Hess S_W(v,v) ≥ m|v|²`. At `β=2.5`: `Ric` of the group factor is a `β`-independent positive constant (SU2 `½g`, SU3 `¾g`); `Hess S_W ∝ β` reaches `≈−1.1β=−2.8` (SU2 L=4 min `−2.81`, 21.7% of 3072 directions negative; SU3 L=3 min `−2.77`, 42%). True Ric ⇒ `min(Ric+Hess S)≈−2.3/−2.0`; even a generous `Ric=+1.5` gives `−1.3`. **So `m<0` for every `m>0`.** The criterion does not merely resemble the Savvidy-sign wall — it *is* `Hess(−D²+2adF)` again (PROOF-SAVVIDY-EXACT-SIGN: `λ_min=−g|B|`, negatives real physics, W4 measures 7.1% gauge-projected transverse negative at the SU(2) crossover after the quotient). BE only holds in `β ≲ Ric/1.1 ≈ 0.5–1.4` — the SZZ strong-coupling window (`β<1/48` in d=4), i.e. exactly where P3/SZZ already gives the gap.

**(b) Coarse space = H1, losslessly (zero slack — the killer).** B3-FUNCTIONAL-BE §1: the coarse generator lives on the flat pushforward of `P` (an orthogonal projection onto a linear subspace of `X=log U`). Take `f` **affine** (`Hess f≡0`, killing BE's *only* possible positive term `‖Hess f‖²`):
```
   Γ₂(v,v) = Hess S_eff(v,v) = G″(v,v).
```
There is **no loss**: `Poincaré(ρ)` on `ν` is *equivalent* to `G″⪰ρ`, hence to `m_eff²(k→0)≥ρ`, hence to `H1` (OTHER-H1-POINCARE). Flatness removes the entire advantage BE could have over the trivial identity `Hess V_eff=G″`. So "prove `Γ₂≥ρ`" and "prove `H1`" are *the same sentence* — the machinery is the wrong half (PROOF-IR-H1 in spirit), not a lever.

**(c) The "proven Haar part" is the wall's own slack term, ~10× short.** The proposal's `Ric_Haar`/`W_geo=E[Φ″]−Var(Φ′)≥0` (Prékopa, PROOF-JACOBIAN-SAVES) is a **genuine theorem** — every mode incl `k→0`, every `L`, value-free. But it is precisely the `W_geo` already banked inside `G″=Δ_dyn+W_geo+X`, and it supplies only the `|F|`-free floor `Φ″⪰N/12` (SU2 0.167, SU3 0.25) against `−H_ff≈−1.12β`. No gauge-invariant norm makes the positive part dominate with `c<1`: spectral ratio `ρ≈3–16` at `β=2.5`, and the only `β` admitting pointwise `c<1` is `β<0.15` (PROOF-SAVVIDY-WGEO-DOMINATION). The extra convexity beyond `N/12` (`g″→∞`) sits only at the log-cut `d→2π` (links→center `z`, `F≈0`) — **anti-correlated** with the Savvidy driver `∝|F|` (B3 §3(iii)). So "Haar proven + action rescued" decomposes the wall into [banked slack] + [the wall]; it does not remove either.

**Verdict on Kill 1: LANDS.** `E[Γ₂]≥ρ` is the Savvidy-sign / `Hess≥λ` wall, renamed. Beat #2 (sign-controlling lower bound) FAIL.

---

## §2. KILL 2 — integrated-vs-pointwise illusion: measure-starvation starves the wrong tail

**Lands — decisively; this is the proposal's own chosen escape hatch, and it is already proved insufficient.**

The proposal's hope is `E[Γ₂]≥0` despite pointwise negatives, via WSTC (the Ward–Stein transport cone bound, ABC-NEW-INEQ). Two independent horns both close it.

**(a) The starvation targets the wrong region.** WSTC proves the **high-field cone** `C_κ={‖F‖²≥κ²}`, `κ>κ_typ`, contributes only
```
   Var_{μ_0}(f·1_{C_κ}) ≤ (cβκ_max)² · e^{−cβ(κ²−κ_typ²)N_□} → 0,   uniformly in L.
```
But the negative variance that decides the sign is **not** high-field. W6-BE §3 and ABC-NEW-INEQ §4: at `β=2.5` the negatives are **high-momentum at typical field** (SU2 `−1.1β` at 21.7% of directions, SU3 42%, **90–97% in the fluctuation/UV block**, W4 7.1% transverse-negative at the SU2 crossover). These sit at `‖F‖≈κ_typ`, **outside `C_κ` for every `κ>κ_typ`**, where WSTC's exponent vanishes. Measure-starvation removes the negative variance where there is no physical instability (high field is already `e^{−β‖F‖²V}`-suppressed) and leaves untouched where the instability actually lives. It is a correctly-signed theorem about the **wrong tail** — ABC-NEW-INEQ's exact verdict: "WSTC bounds the wrong tail; the residual is exactly H1."

**(b) The soft `k→0` marginal stays negative absent conditioning.** W6-BE §3: "At `k→0`, where blocked/induced mass would enter, the restoring `p²→0` and the marginal proxy `G″ ≤ E[PᵀHP]` soft value `≈+0.043−0.208=−0.165<0`." The coarse Savvidy mode `PᵀHP` min `≈−0.3` (MARGINAL) survives blocking at full strength because blocking is a low-pass filter and the instability is *at* `k→0`; `Var≥0` only subtracts (`G″≤E[PᵀHP]`). The measured `+0.51`/`+0.27` positive `G″` is carried by the **conditioning** `E[H|0]=+0.22` and `W_full≈+0.54` — a *measure* fact about `Var(G|0)`, which is H1 (PROOF-JACOBIAN-SAVES, OTHER-H1-POINCARE).

**(c) Even granting the integration, it reduces to H1.** On the coarse mode the exact Ward identity `G″·Var(a)=1` (ABC-NEW-INEQ §1) ties the curvature to the variance: `E_ν[Γ₂]=G″=E[H|0]−Var(G|0)=m_eff²`, so `E[Γ₂]≥ρ ⟺ Var_μ(a_{k→0})≤1/ρ = H1`. There is no "integrated" statement distinct from the gap. (This is the task's second horn — "`E[Γ₂]≥ρ` still requires the gap" — and it is exact, not merely plausible.)

**(d) A correction to the proposal's premise.** The task's phrasing "the softest unstable mode is at LOW field / extremal-antiselfdual configs" needs the corpus's correction (PROOF-SAVVIDY-EXACT-SIGN §3): extremal (BPS, SD or ASD) backgrounds are **PSD** (`λ_min=0`, kernel=moduli) — they are the *safe* configs. Negativity needs **non-extremality**, maximal at `Q=0` (the equal SD/ASD Savvidy mixture). So the picture is: the single deep `k_∥=0` Savvidy rung `−g|B|` sits on a **large-field** non-BPS saddle (which WSTC *does* starve, but it is measure-suppressed), while the **populated** negative variance is typical-field/high-momentum (which WSTC does **not** starve). Either way the rescue is on the wrong object. **The "starvation removes them" claim is false.**

**Verdict on Kill 2: LANDS.** The negative modes are not measure-starved in the relevant direction; WSTC suppresses the wrong region.

---

## §3. KILL 3 — free-field leak / vacuity: it certifies a spectral gap, not `f₁²>0`, and the "Haar" part is a compactness artifact

**Lands.**

**(a) Spectral gap ≠ atom.** Bakry–Émery certifies the Poincaré / LSI constant of the transfer operator — a *disorder* statement `Var_μ(φ)≤(1/λ)ℰ_μ(φ)`. That is exactly `H1`, which the campaign has established is **necessary but not sufficient** for the Clay object: `f₁²>0` does not imply `m>0` (MS-BLIND §3: "a massless excitation carries `f₁²>0` at a pole `s=0`; the residue is the intensity; the gap is its lower support edge"), and a variance/Poincaré bound is **blind to atoms** — the whole class of positivity tools (RP/LSI/FKG/BL/Prékopa) "only certifies disorder and is blind to atoms" (SIGNED-FOURIER-BLIND §0). So the machinery does not "natively give a gap lower bound" for the Clay object; it gives `H1`, which *is* the unproven step. Proving `Γ₂≥ρ` = proving `H1` = proving the variance bound; it does not prove `f₁²>0` (which additionally needs `I2`, isolation — PROOF-P3 §5 corrects any claim that grading yields `I2`).

**(b) The free-field acid test leaks on the "massive" side.** On the free **massless** continuum field, `Φ=0`, `Hess S=0`, zero mode flat ⇒ `Γ₂=0`, returns `0` — the acid test "passes" trivially (`0≥0`, INVENT-NEW-MONOTONE §4(i), B3 §5). But the same machinery returns `m²>0` on the free **massive** field. That is exactly the failure the task names: it **cannot distinguish** "free massive gap" from "YM glueball atom" — it reports whatever dynamical mass is already in the measure and does not manufacture the interacting atom.

**(c) The "proven-nonnegative Haar part" is the sink.** The proposal's pride — `Ric_Haar≥0`, `W_geo≥0` — is on the lattice the floor `N/12` that is **present even at `β=0`** (`0.248`, B3 §5 caveat). `Ric` of compact SU(N) is a `β`-independent positive constant (`½g`/`¾g`, W6 §1): it is the **curvature of the group manifold**, i.e. a `k²→0` finite-volume / compactness artifact, absent for the continuum free field. Attributing the gap lower bound to the Haar part therefore certifies a **kinematic** mass that vanishes in the continuum — a free-field leak in the compactness direction. To avoid it one must attribute the sign to the **conditioning** `E[H|0]=+0.22` and `W_full≈+0.54` — i.e. to `Var(G|0)`, i.e. to H1 (PROOF-JACOBIAN-SAVES, OTHER-SAVVIDY-SIGN).

**Verdict on Kill 3: LANDS.** BE gives the spectral/Poincaré gap (= H1), blind to the atom; the "proven Haar part" is either a compactness artifact (leak) or a restatement of the conditioning (`Var(G|0)` = H1). Beat #4 FAIL (with leak).

---

## §4. KILL 4 — circularity: `ρ~Λ²` is the scale/mass you are deriving

**Lands.**

To get `Γ₂≥ρ` at the physical point you need `ρ>0` with `ρ~m_eff²~Λ²`, and the exact missing inequality is (◆) of PROOF-JACOBIAN-SAVES:
```
   Var_μ(⟨v,∇S⟩|a=0) − E_μ[H_vv|a=0] < W_full,   uniform (L,U′)
   ⟺ G″⪰ρ ⟺ Var(a_{k→0})≤1/ρ ⟺ m_eff²(k→0)≥ρ.
```
- **Route B (BE on `ν`):** `Hess V_eff⪰ρ` *is* H1 verbatim (OTHER-H1-POINCARE route B) — you assume the bound to prove it.
- **Route A (push full LSI through `P`):** `Var_ν(φ)=Var_μ(φ∘P)≤(1/λ_μ)ℰ_μ(φ∘P)` needs `λ_μ` = the full `β=2.5` gap (proved only `β<1/48`) — circular. W6 §3: "the induced mass would be the gap itself — circular." The corresponding optimal-control/HJB precedent is the same: "a lower bound needs the inequality uniform over controls = LSI/H1; the optimizer is `ψ₁`" (REDTEAM-FORGE, PROOF-OPTIMAL-CONTROL-HJB).
- **Reweighting** (Holley–Stroock bridge from pure-Haar to `μ`) inflates the LSI constant by `e^{Θ(N_□)}` (OTHER-SAVVIDY-SIGN §2) — vacuous, and the only non-vacuous regime is `β<1/48` where SZZ already holds.

The bound's constant `ρ` is the gap; supplying `ρ~Λ²` from outside is supplying the answer. **Kill 4: LANDS.** Beat #3 FAIL.

---

## §5. Five-point verdict (honest)

| # | Constraint | Verdict | Exact failure |
|---|---|---|---|
| 1 | Nonperturbative | **PASS** | BE/Prékopa/WSTC are measure/geometry-level; WSTC and `W_geo≥0` are genuine nonperturbative theorems. |
| 2 | Sign-controlling **lower** bound `m_eff²(k→0)≥ρ` | **FAIL** | `Γ₂=G″` losslessly on the coarse flat space (B3 §1); full-space `Ric+Hess S≈−2.3` (W6 §1); Haar part `N/12` ~10× short and anti-localized (WGEO-DOMINATION); WSTC starves the wrong tail (ABC-NEW-INEQ §4). |
| 3 | Non-circular | **FAIL** | `Γ₂≥ρ` ⟺ `G″⪰ρ` ⟺ H1 verbatim (route B); route A needs the full `β=2.5` LSI (route A circular); W6 "induced mass = gap itself." |
| 4 | Non-SUSY + **fails free field** | **FAIL (leak)** | Continuum free massless ⇒ 0 (fine) but machinery returns free mass on free massive (cannot distinguish atom); the "proven Haar part" is the compactness artifact `N/12` present at `β=0` — a kinematic mass, not the gap. |
| 5 | Uniform `≥ρ`, reach `β≈2.5` | **FAIL** | Exactly H1; only `≥0` uncontrolled at 2.5, genuine `≥N/12` only at `β=0`, `λ=C_blk^{−2}λ_μ` only at `β<1/48`. |

**Survival score: 1/5** (only nonperturbativity). Identical profile to SIGNED-FOURIER-BLIND (1/5) and the INVENT-NEW-MONOTONE / W6-BE / B3-FUNCTIONAL-BE precedents.

---

## §6. Why this is a re-tread, not a new species (the decisive prior art)

The species is the **union of three already-dead probes**, none in a new class:
- **W6-BE-at-2.5** — full-space `CD(m,∞)`/`Ric+Hess S_W≥m`: FALSE, `min≈−2.0/−2.3`.
- **B3-FUNCTIONAL-BE** — coarse/conditioned `Γ₂`: *exactly* H1, lossless, zero slack; dies at the `k→0` cone mode.
- **INVENT-NEW-MONOTONE** — the Perelman `Γ₂`-refinement as the (only) gap-certifying monotone: "`Γ₂=Ric_A+Hess_A S_W` is Savvidy-indefinite, `min≈−2.0`"; "the `Γ₂`-refinement — the *only* gap-certifying monotone — is sign-indefinite exactly where the gap must be proved." REDTEAM-FORGE already lists INVENT-NEW-MONOTONE as breaks {3,5}: "the monotone that holds is relative entropy to `μ` — circular; the gap-certifying `Γ₂` is sign-indefinite at 2.5."

The proposal's one addition — "Haar part proven + Savvidy part rescued by WSTC" — is not a new mechanism: the "Haar part proven" is the banked `W_geo≥0` (already the slack term of `G″`), and "rescued by WSTC" is *provably the wrong tail* (ABC-NEW-INEQ). So the decomposition is exactly `G″ = W_geo + (E[H]−Var)`: [proved slack] + [H1]. The optimal-transport/Riemannian framing adds no coercivity the variance identity `G″Var(a)=1` lacks — on flat space it is literally the same quadratic form.

---

## §7. What, if anything, survives (the honest residue)

- **Genuine, already-banked, not new:** `W_geo=E[Φ″]−Var(Φ′)≥0` (Prékopa) — a real theorem on every mode incl `k→0`, every `L` (PROOF-JACOBIAN-SAVES). Also the shallow-cell closure `Γ₂≥0` on `{2g|F|≤N/12}` (B3 §3(i), Brascamp–Lieb) and the high-field WSTC bound — all real, all already proved, none touching the open cell.
- **A faithful dialect, not a closure:** `Γ₂` on the flat coarse space *is* `G″` (B3 §1). This is a useful restatement of H1 in optimal-transport language — it lets one cite `HWI`/`CD`/`LSI` machinery — but it is **lossless**, so it cannot weaken H1; "H1 ⇒ gap, H1 unreduced" (REDTEAM-FORGE cross-cutting).
- **The exact missing OBJECT is one conditioned-variance bound** (◆): `Var_μ(⟨v,∇S⟩|a=0) ≤ E[V″|a=0]−ρ` at typical field — the same single "variance upper bound" species the skeleton names at G2.5/H1 and G1.2/L3 (docs/PROOF-SKELETON.md §0, §3). The transport framing does not supply it; it re-expresses the demand for it.

---

## §8. Sharpest line

**The Bakry–Émery/optimal-transport species dies on the one fact the campaign already banked: on the flat coarse space where the gap lives, `Γ₂` is lossless and equals `G″` exactly — so "`Γ₂≥ρ`" *is* `H1` with zero slack, not a lower bound on it (B3-FUNCTIONAL-BE §1). Its "proven Haar part" is the banked slack `W_geo≥0` (`N/12`, ~10× short of the Savvidy `−1.1β`, and anti-localized at the log-cut where `F≈0`); its "rescue by measure-starvation" (WSTC) provably starves the high-field cone, where the negative weight is measure-suppressed anyway, and leaves untouched the sign-deciding negative variance at typical field and high momentum (W6-BE §3, ABC-NEW-INEQ §4) — the wrong tail. Full-space BE is the Savvidy wall verbatim (`Ric+Hess S≈−2.3<0`, W6 §1); the coarse BE is H1 verbatim; the Perelman `Γ₂`-monotone is sign-indefinite exactly where the gap must be proved (INVENT-NEW-MONOTONE). The optimal-transport framing adds no coercivity the variance identity lacks — it is a faithful rename of `G″=W_geo+(E[H]−Var)`, [proved slack] + [H1], not a lowering of the wall.**

**Exact stuck step.** `Var_μ(⟨v,∇S⟩|a=0) ≤ E_μ[V″|a=0]−ρ` at typical field, uniform in `(L,a)`, with `ρ>0` — (◆) of PROOF-JACOBIAN-SAVES = `Var_μ(a_{k→0})≤1/ρ` = H1. Producing it *is* the missing "variance upper bound" species; the Γ₂/OT calculus does not produce it, it *is* the demand for it.

---

*Refs.* Internal (read/cited as ground, not invented): W6-BE-at-2.5 (`CD(m,∞)` attempted, `Ric+Hess S≈−2.3/−2.0`, `−1.12β`, `m<0`; false); B3-FUNCTIONAL-BE (§1 `Γ₂=G″` lossless = H1; §3 anti-localization; §5 `N/12` compactness caveat); OTHER-H1-POINCARE (route A circular / route B = H1; blocking is low-pass, `PᵀHP≈−0.3`); PROOF-COARSE-POINCARE (exact coarse Dirichlet form, `G″` Savvidy-indefinite, Lemma SC at `β<1/48`); PROOF-JACOBIAN-SAVES (`G″=Δ_dyn+W_geo+X`, `W_geo≥0` Prékopa, (◆)); PROOF-SAVVIDY-WGEO-DOMINATION (no GI norm, `c<1` false, `ρ≈3–16`, log-cut anti-correlation); PROOF-SAVVIDY-EXACT-SIGN (`−D²+2adF`, `λ_min=−g|B|`, negativity non-BPS, max excess at `Q=0`); W4-blocked-positivity (7.1% gauge-projected transverse negative at SU2 crossover — real physics); ABC-NEW-INEQ (WSTC starves high-field cone, wrong tail; negatives typical-field/high-momentum; Ward `G″Var(a)=1`); OTHER-SAVVIDY-SIGN (parity `Q=0` balanced; Holley–Stroock `e^{Θ(N_□)}`); INVENT-NEW-MONOTONE (Perelman `Γ₂` monotone sign-indefinite at 2.5; free `λ_min→0`); REDTEAM-FORGE (INVENT-NEW-MONOTONE breaks {3,5}; cross-cutting finite-regulator acid-test flaw; "H1 ⇒ gap, H1 unreduced"); SIGNED-FOURIER-BLIND (§0 positivity tools blind to atoms; `f₁²` vs mass); MS-BLIND (§3 residue ≠ gap; massless pole carries `f₁²>0`); docs/PROOF-SKELETON.md (§0/§3: single species = variance upper bound at G2.5/H1 and G1.2/L3); docs/YANG-MILLS-MASTER.md lines 117,145. External (pre-vetted in-repo): Bakry–Émery 1985 (structure `Γ₂`); Otto 2001; Sturm 2006; Lott–Villani 2009; von Renesse–Sturm 2005 (CD ⇒ LSI ⇒ gap); Prékopa 1971; Brascamp–Lieb 1976; Holley–Stroock 1987; Savvidy, PLB **71** (1977) 133; Nielsen–Olesen, NPB **144** (1978) 376; Shen–Zhu–Zhu, CMP **400** (2023) 805, arXiv:2204.12737.
