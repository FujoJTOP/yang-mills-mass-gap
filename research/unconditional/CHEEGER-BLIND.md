# CHEEGER-BLIND — independent blind red-team on the "Cheeger isoperimetric spectral-gap inequality" species for the 4D pure SU(N) mass gap

> 2026-11-23 · **判词：SPECIES DEAD — and ALREADY-BANKED. This is not a new species; it is the Cheeger/conductance face of the SHAPE-B / C-BBS / capacity family, blinded to death earlier as `B3-ISOPERIMETRIC`, `LOC-SUFFICIENCY-CHECK`, `RCD-SYNTHETIC-RICCI`, `POTENTIAL-THEORY-CAPACITY`, `C-BBS-FORMULATE/-ARRHENIUS/-SUFFICIENCY`.** The proposal claims the Cheeger inequality `λ₁ ≥ h²/4` escapes the Bakry–Émery/`Γ₂` wall because *isoperimetry needs no curvature hypothesis* (true — Cheeger's inequality is hypothesis-free), so it "stays positive under negative curvature where `Γ₂<0`." That single observation is correct and is the species' only genuine edge. It is also a non-lever: the one object Cheeger demands is `h>0` **uniformly in `L`**, and for this measure `h>0` uniform **is** the coarse-mode Poincaré `λ₁>0` uniform **is** `Var_ν(a_{k→0})≤1/λ` **is H1** — with the diagonal `h=√(λ₁)` on the coarse line and `λ₁=G″=1/χ` on the marginal, so "prove the Cheeger constant" and "prove the gap" are *the same sentence*. Kill 1 lands (on its correct disjunct), Kills 2/3/4 land. **Survival score 2/5** (only nonperturbativity + the non-vacuous free-field acid test, both shared with the banked capacity species). Exact stuck step: a *fixed-scale, `L`-uniform* conductance floor `Φ_ν≥c` on the coarse marginal = coarse-mode localization = `m>0` = H1.

---

## §0. The tool, stated precisely (kill the real thing, not a strawman)

On the coarse zero mode `a=⟨v,X⟩` (block-averaged gauge amplitude), marginal `ν(da)=p(a)da`, `p∝e^{−V̄}`, reversible Dirichlet form `ℰ_ν(f)=∫|f′|²dν`, generator `Lf=f″−V̄′f′`, Poincaré constant `λ₁(ν):=inf_{Var f>0}ℰ_ν(f)/Var_ν(f)`. The proposal:

```
Cheeger:   λ₁ ≥ h²/4,   h(ν) = inf_A  ∂_ν A / min(ν(A),ν(A^c))   [Cheeger 1970, hypothesis-free]
target:    λ₁>0  ⇒  Var_ν(a) ≤ ℰ_ν(a)/λ₁ = 1/λ₁   (ℰ_ν(a)=1)  ⇒  G″(k→0) ≥ λ₁  (Ward G″·Var=1)
claim:     isoperimetry stays positive under negative curvature (no Γ₂=Ric+Hess S needed),
           so it may cross the Savvidy crossover where Bakry–Émery (Γ₂) cannot.
```
Everything below tests five constraints: (1) nonperturbative; (2) a sign-controlling **lower** bound `m_eff²(k→0)≥λ`; (3) non-circular; (4) fails the free field (`f₁²=0`,`G″=0`); (5) `(L,a)`-uniform, reaches `β≈2.5`.

**Why this is a re-tread (Step 0, the prior version).** The Cheeger/isoperimetric route is already in-repo as `B3-ISOPERIMETRIC` (title: "building the coarse Poincaré constant ρ **directly** via Cheeger/isoperimetry on ν"), and its refinements are `LOC-SUFFICIENCY-CHECK` (1D Cheeger constant = global bottleneck), `RCD-SYNTHETIC-RICCI` (Buser `ρ(K,D)∼1/D²→0`), `POTENTIAL-THEORY-CAPACITY` (`λ₁≥Φ_eq²/4`, capacity floor = H1), `C-BBS-FORMULATE/-ARRHENIUS/-SUFFICIENCY` (conductance floor = marginal Poincaré = H1). The 9th species is the union of these, re-dressed. The attack below is on *that* prior version, with the corpus's own numbers.

---

## §1. KILL 1 — Cheeger **is** H1 renamed (`h>0` uniform ⟺ `λ₁>0` uniform)

**Lands — on the correct disjunct. (The other stated mechanism is mis-aimed; I report which.)**

The task offers two mechanisms. Be honest about each:

**(a) "the isoperimetric profile at the crossover collapses; concentration of `e^{−β‖F‖²V}` on near-minimal `|F|` drives the boundary fraction over negative-mode directions to 0" — DOES NOT LAND AS STATED.** The corpus refutes this picture on two counts. First, the sign-deciding negative variance is **not** measure-starved: `W6-BE`/`ABC-NEW-INEQ` (§4) put it at **typical field, high momentum** (SU2 `−1.1β` at 21.7% of directions, SU3 42%, 90–97% in the fluctuation/UV block; W4 7.1% gauge-projected transverse-negative at the SU2 crossover), *inside* the thermal shell where the mass lives, not on a starved cone. Second, the coarse isoperimetric profile does **not** collapse to 0 — the measured thermal-shell barrier is `H∼|V̄″|λ⁻¹≈0.3/λ=O(1)` at the target scale (`POTENTIAL-THEORY-CAPACITY` §3, `C-BBS-ARRHENIUS` §3). The profile is `O(1)`, not `→0`. A collapsing profile would actually be the **danger**, not the escape.

**(b) "Cheeger's constant for this measure is numerically = the gap" — LANDS, decisively.** On the 1D coarse line the Cheeger constant is *by definition* the global bottleneck of `ν` (`LOC-SUFFICIENCY-CHECK` §1–2: one-dimensional Cheeger constant = minimum bottleneck over the whole axis, taken where the isoperimetric profile is most expensive). Via Cheeger–Buser `h≍√λ₁`, and on the marginal `λ₁=1/χ=G″` (Gaussian-exact Ward `G″·Var=1`; `D3`, `OTHER-H1-POINCARE`), the number is literally the gap: measured `G″=+0.51` SU3, `E[H|0]=+0.22`, `W_full≈+0.54` (`W7` Part B; `OTHER-H1-POINCARE` §accounting). So **"prove `h>0`" is "prove `λ₁>0`" is "prove `G″≥λ`" is H1 — zero slack.** The corpus states this exactly (`B3-ISOPERIMETRIC` sharpest line):

```
ρ ≥ h(ν)²/4,  h(ν)=inf_B ∂_ν B / min(ν(B),ν(B^c)),  B fiber-saturated;
∃ρ>0 uniform in L  ⟺  h(ν) ≥ 2√ρ  ⟺  transverse bottleneck of ν bounded below
                   ⟺  Var_ν(a_{k→0}) ≤ 1/λ  ⟺  m_eff²(k→0) ≥ λ   (H1).
```
And the two-coordinate obstruction (`B3-ISOPERIMETRIC` §"new no-go") kills the *only* thing that could make Cheeger easier than `λ₁`: the fiber's own LSI/isoperimetry is **silent** on the `ν`-Cheeger competitors `P̄⁻¹(B)` (fiber-saturated sets carry identically zero fiber-tangent boundary), so no fiber tool can supply `h(ν)>0`. The marginal `ν` and the conditionals `{μ(·|a)}` are independent coordinates of `μ=ν·μ(·|a)`.

**Verdict on Kill 1: LANDS** (disjunct b; disjunct a mis-stated — corrected here). "Cheeger isoperimetry" is the H1 wall wearing a hypothesis-free inequality, not a way past it. Beat #2 FAIL.

---

## §2. KILL 2 — Cheeger→Poincaré does **not** close; the local-Poincaré input is the wall

**Lands.**

Cheeger `λ₁ ≥ h²/4` is the direction that *would* help (a lower bound on `λ₁` from a conductance floor, no curvature needed). But it is only usable if `h>0` is known; and the object it gives is a **one-sided** bound. The reverse (Buser `λ₁ ≤ C(h+√h)`, giving the two-sided `h≍√λ₁` that would make the equivalence robust) needs a **lower Ricci bound** — which is exactly the Savvidy-falsified `Γ₂` (`Ric+Hess S≈−2.3`, `W6-BE`). So the "isoperimetry escapes curvature" hope is a fork with both horns dead:

- **Use only the forward Cheeger inequality:** then you must *supply* `h>0`, and the constructive way to lower-bound `h` (a conductance floor / equilibrium-capacity / local-Poincaré-on-balls) is a **local Poincaré inequality**, which on the loop/gauge space is precisely the missing coarse-mode Poincaré (`C-BBS-FORMULATE` §4 (R4); `C-BBS-SUFFICIENCY` §1). `λ₁ ≥ Φ_eq²/4 ⟹ χ ≤ 4/Φ_eq²` — verbatim H1.
- **Try the two-sided Cheeger–Buser:** it demands `Ric ≥ −K`; the uniform `L→∞` constant collapses regardless — `ρ(K_*,D) ∼ c(K_*)/D² → 0` as `D∼L²→∞` (`RCD-SYNTHETIC-RICCI` §4), a finite-`V` kinematic artifact, not `m_eff²(k→0)≥λ`.

The one genuine bonus — Cheeger's inequality is **hypothesis-free**, so it does not need `Γ₂≥0` and is therefore not *literally* the `Γ₂`-sign wall — is real but non-load-bearing: it relocates the unknown from `λ₁` (unproven) to `h` (equally unproven, and equal to it). You cannot lower-bound `h` without the local Poincaré that is H1. **Verdict on Kill 2: LANDS.** Beat #5 FAIL; the exact missing input is a local (ball) Poincaré / conductance floor uniform in `(L,a,scale)`.

---

## §3. KILL 3 — free-field / vacuity: certifies the **gap** (H1), not the **atom** `f₁²`; no `>0` leak, but atom-blind

**Lands (with one honest credit).**

**(a) The acid test is genuine — non-vacuous, the species' best feature, shared with capacity.** Free massless: zero mode flat/non-compact, `V̄≡0`, `Cap=0`, `h=0`, `λ₁=0`, `Var(a)=∞` (`POTENTIAL-THEORY-CAPACITY` §4; `C-BBS-SUFFICIENCY` §3; `RCD-SYNTHETIC-RICCI` §5). So Cheeger returns **no** gap on the free field — a real discriminator (a flat-direction moment/conductance), stricter than `B3`'s fiber-LSI which was blind to the free transverse mode. It does **not leak `>0`** on the free field. Credit where due: this is the cleanest non-vacuous acid test in the whole positivity family. But see (c) — the acid-test pass can be bought with a *definitional swap*.

**(b) What it certifies is `H1` (necessary), not the Clay atom `f₁²>0`.** Cheeger bounds a Poincaré / second-moment datum — a *disorder* statement (`Var_ν(φ)≤λ⁻¹ℰ_ν(φ)`), which is exactly the necessary H1. The atom `f₁²>0` lives at the **support edge** `inf supp ν_O(·)`, a two-level (`λ₁−λ₂`) datum that every bottom/Dirichlet/Poincaré/Cheeger functional is **blind to** (`C-BBS-SUFFICIENCY` §3, §"Scorecard" #5; `PROOF-EDGE-E2`; `G2G-TRANSFER-SPECTRAL` §4). The audit counterexample is decisive: `ν_O=wδ_ε+(1−w)δ_S`, `w/ε≤1/λ` ⇒ finite Cheeger/capacity constant (`λ₁>0`) **while** `inf supp=ε→0` (no gap). A Cheeger floor, uniform or not, cannot rule this out — it sees `χ`, not the thin `δ_ε` mass. Cheeger one rung too low, exactly as capacity: **necessary, not sufficient.**

**(c) The circular/sufficiency leak.** To convert `λ₁>0` into the physical mass gap `Δ=inf supp` one must already assume the residue/isolation structure (`f₁²>0` with isolation `I2 = H*`) — that is the circular input. And the *only* non-vacuous version uses the **global** Cheeger constant `h` (which is H1); the "local single-cell oscillation" version is the definition-swap already nailed in `LOC-SUFFICIENCY-CHECK` §3 (`H₀:=−logΦ_ν` global vs local cell) and `C-BBS-ARRHENIUS` §4 — the local reading self-destructs (`T(λ)≤λ/7`, contracts to `λ=0`, `C-BBS-ARRHENIUS` §4). So either you certify H1 (insufficient), or you take the global `h` (= H1, unproven), or you leak through the definition swap. No fourth option.

**Verdict on Kill 3: LANDS.** BE's exact failure (§3a of `TRANSPORT-BLIND`) repeats: the species certifies the spectral/Poincaré gap, blind to the atom; the "cheapest acid test" is genuine but stops at H1. Beat #4 PASS, but it buys only a necessary condition.

---

## §4. KILL 4 — circularity: to bound `h` below you must know the global concentration of `μ` (the gap's own statement); `h∼√Λ`

**Lands.**

To get `h>0` at the physical point you must lower-bound the conductance of **every** cut of `ν` — i.e. know the global confining concentration of the measure. On compact finite volume `h>0` is trivial (`ν` compactly supported ⇒ `Var<∞`, `Cap>0` for *any* measure, `POTENTIAL-THEORY-CAPACITY` §2; `C-BBS-ARRHENIUS` §1); the entire content is `L`-uniformity, and the `L`-uniform floor is **gap-strength**. Concretely (`POTENTIAL-THEORY-CAPACITY` §exact criterion):

```
h>0 uniform  ⟺  V̄ coerces ½λa² at the tail  ⟺  ∫(1+a²)e^{−V̄}da<∞ uniformly  ⟺  Var(a)≤1/λ  (H1),
```
with the constant `h ∼ √λ₁ ∼ √Λ` — the mass scale you are deriving. The Savvidy dip on the thermal shell is precisely the barrier whose finiteness controls `h`: `λ₁ ≳ e^{−2H}`, `H∼0.3/λ` — `O(1)` at the target (`C-BBS-ARRHENIUS` §3), so no exponential catastrophe, but proving `H` is `L`-uniform is, verbatim, coarse-mode localization (`C-BBS-ARRHENIUS` §5: the exact residual is "equivalent to localization of the coarse mode (= `m>0` In capacity costume)"). Supplying `h>0` from outside is supplying the answer. **Verdict on Kill 4: LANDS.** Beat #3 FAIL, `ρ~Λ` circular.

---

## §5. Five-point verdict (honest)

| # | Constraint | Verdict | Exact failure |
|---|---|---|---|
| 1 | Nonperturbative | **PASS** | Cheeger/coarea/capacity are measure-level, all `β`, no loops. |
| 2 | Sign-controlling **lower** bound `m_eff²(k→0)≥λ` | **FAIL** | `h>0` uniform ⟺ `λ₁>0` (1D Cheeger constant = global bottleneck, `LOC-SUFFICIENCY-CHECK` §1); trivially `>0` on compact `L`, content is `L`-uniformity = H1. Hypothesis-free inequality (real) but non-load-bearing. |
| 3 | Non-circular | **FAIL** | `h>0 ⟺` tail-coercivity `½λa² ⟺ Var(a)≤1/λ` = H1 (`POTENTIAL-THEORY-CAPACITY` §2); Buser two-sided needs `Ric≥−K` (`RCD` §4, `ρ∼1/D²→0`); two-coordinate obstruction blocks fiber supply (`B3-ISOPERIMETRIC`). |
| 4 | Non-SUSY + **fails free field** | **PASS** | Flat free zero mode ⇒ `Cap=0`,`h=0`,`λ₁=0`,`Var=∞` — cleanest non-vacuous acid test; **but** certifies `H1`, atom-blind (`C-BBS-SUFFICIENCY` §3 `wδ_ε` counterexample). |
| 5 | Uniform `≥ρ` reach `β≈2.5` | **FAIL** | Only H1; compactness floor dies `ρ(K,D)∼1/D²→0` (`RCD` §4); genuine `≥N/12` (0.167 SU2/0.25 SU3) only at `β=0`, `λ=C_blk⁻²λ_μ` only at `β<1/48` (`Lemma SC`, SZZ). |

**Survival score: 2/5** — identical profile to `POTENTIAL-THEORY-CAPACITY` (2/5) and `C-BBS-SUFFICIENCY` (necessary, not sufficient). Passes only nonperturbativity (1) and the non-vacuous free-field acid test (4); both are already banked in the capacity species. Fails the Clay object on (2),(3),(5).

---

## §6. Why this is a re-tread, not a 9th species (decisive prior art)

The proposal is the **Cheeger face of the already-blinded capacity/Dirichlet family**, none in a new class:
- **`B3-ISOPERIMETRIC`** — title literally "Cheeger/isoperimetry on ν"; the two-coordinate obstruction (fiber LSI silent on `P̄⁻¹(B)`), the `h(ν)>0 ⟺ H1` equivalence, verdict "No constructive ρ exists."
- **`LOC-SUFFICIENCY-CHECK`** — 1D Cheeger constant = **global** bottleneck; single-cell `O(1)` oscillation ⇒ global floor is **NO-GO** (free-field + double-well `a=±a₁` counterexamples); Cheeger/Poincaré algebra = true theorem, but the local→global upgrade = H1.
- **`RCD-SYNTHETIC-RICCI`** — exactly the Cheeger–Buser/`RCD` framing: on flat `A`, synthetic `CD(K,∞) ⟺ Hess V_eff⪰K`; the only positive Poincaré is compactness `ρ(K,D)∼1/D²→0`.
- **`POTENTIAL-THEORY-CAPACITY`** — `λ₁≥Φ_eq²/4` (Cheeger/Sinclair); exact tail-coercivity criterion; capacity floor = H1; 2/5.
- **`C-BBS-FORMULATE` / `-ARRHENIUS` / `-SUFFICIENCY`** — the Cheeger **conductance floor** built explicitly (`λ₁≥Φ_ν²/4`, "no indefinite operator is signed"); the floor is the marginal Poincaré = `1/χ=G″` = H1; necessary, not sufficient.

The proposal's one addition — "isoperimetry needs no curvature, so it escapes `Γ₂<0`" — is **true and already noted** (`C-BBS-FORMULATE` §1: "the quadratic form is `⩾0` for *any* `V̄`; no indefinite operator is signed"; `POTENTIAL-THEORY-CAPACITY` §1: "trades 'sign an indefinite operator' for 'min a positive energy'"). It is a *cosmetic* gain (a sign is obtained for free) that leaves the *quantitative* lock (`h>0`/`Φ≥c` uniform `L`) untouched. Same lock, cleaner-looking key.

---

## §7. What, if anything, survives (the honest residue)

- **Genuine, already-banked, not new:** the hypothesis-free sign (`λ₁≥h²/4≥0`, `λ₁≥Φ_eq²/4≥0` — positivity for free, no `Γ₂`), and the **non-vacuous free-field acid test** (flat free mode ⇒ zero capacity/Cheeger constant — a genuine geometric discriminator absent from `B3`'s fiber route). Both are real; both are `POTENTIAL-THEORY-CAPACITY`/`C-BBS` results.
- **A faithful dialect, not a closure:** Cheeger on the coarse line is `h≍√λ₁=√(G″)` — a restatement of H1 in isoperimetric language. It lets one cite Cheeger/Buser/`RCD` machinery; it is **lossless**, so it cannot weaken H1 ("H1 ⇒ gap, H1 unreduced", `REDTEAM-FORGE` cross-cutting).
- **The exact missing OBJECT** is one fixed-scale, `L`-uniform conductance floor / coarse-mode localization bound: `Φ_ν≥c` uniform in `(L,a,scale,k→0)` — the same "variance upper bound" species the skeleton names at G2.5/H1 and G1.2/L3 (`docs/PROOF-SKELETON.md` §0,§3; `docs/YANG-MILLS-MASTER.md` lines 117,145). Cheeger does not supply it; it re-expresses the demand for it.

---

## §8. Sharpest line

**The Cheeger species dies on the one fact its own family already banked: Cheeger's inequality `λ₁≥h²/4` is genuinely hypothesis-free — so it does *not* pass through the Savvidy-falsified `Γ₂`, that part is real — but the object it needs, `h>0` **uniformly in `L`**, is exactly the coarse-mode Poincaré `λ₁>0` (`h≍√λ₁`, with `λ₁=1/χ=G″` the H1 number on the marginal; `B3-ISOPERIMETRIC`, `LOC-SUFFICIENCY-CHECK` §1). On compact finite volume `h>0` is vacuous; the content is `L`-uniformity, and the uniform floor `Φ_ν≥c` is H1 in isoperimetric costume (`λ₁≥Φ²/4 ⟹ χ≤4/Φ²`, `C-BBS-SUFFICIENCY` §1). The two-sided Cheeger–Buser that would make it robust needs `Ric≥−K` — the falsified Savvidy curvature — and its compactness constant dies as `ρ(K,D)∼1/D²→0` (`RCD-SYNTHETIC-RICCI` §4). The free-field acid test is clean (flat mode ⇒ `h=0`), but it certifies a Poincaré/second-moment (H1), blind to the atom `f₁²` (support edge, invisible to any bottom/conductance functional — `wδ_ε` counterexample). Isoperimetry does not cross the wall; it renames `λ₁` as `h` and stops at the same bottleneck.**

**Exact stuck step.** A fixed-scale, `(L,a)`-uniform conductance floor `Φ_ν ≥ c>0` on the coarse marginal — equivalent to `h(ν)≥c` ⟺ `Var_ν(a_{k→0})≤1/c` ⟺ `m_eff²(k→0)≥c` ⟺ coarse-mode localization (`C-BBS-ARRHENIUS` §5) = H1.

---

*Refs.* Internal (read/cited as ground, not invented): `B3-ISOPERIMETRIC` (Cheeger/isoperimetry on ν; two-coordinate obstruction; `h(ν)>0⟺H1`; "No constructive ρ exists"); `LOC-SUFFICIENCY-CHECK` (1D Cheeger constant = global bottleneck; single-cell→global NO-GO; free-field + double-well counterexamples); `RCD-SYNTHETIC-RICCI` (synthetic `CD⟺Hess` on flat `A`; Buser `ρ(K,D)∼1/D²→0`); `POTENTIAL-THEORY-CAPACITY` (`λ₁≥Φ_eq²/4`; exact tail-coercivity criterion; capacity floor = H1; 2/5); `C-BBS-FORMULATE` (Cheeger conductance floor; "no indefinite operator is signed"; (R1)-(R4)); `C-BBS-ARRHENIUS` (thermal-shell `H∼0.3/λ=O(1)`; local reading self-destructs `T(λ)≤λ/7`); `C-BBS-SUFFICIENCY` (`λ₁≥Φ²/4 ⟹ χ≤4/Φ²` = H1; `wδ_ε+(1−w)δ_S` edge≠support counterexample; necessary not sufficient); `OTHER-H1-POINCARE` (coarse form exact; `G″=1/χ`; measured `+0.51`, `E[H|0]=+0.22`); `W6-BE` (`Ric+Hess S≈−2.3/−2.0`, `−1.1β`; negatives typical-field/high-momentum); `ABC-NEW-INEQ` (WSTC starves wrong tail; Ward `G″Var(a)=1`); `MARGINAL` (`PᵀHP` min −0.3); `INVENT-VARIANCE-INEQ` (`D/S≈2%`); `PROOF-SAVVIDY-WGEO-DOMINATION` (`W_geo≥0`, `N/12`, log-cut anti-correlation); `PROOF-COARSE-POINCARE` (Lemma SC `β<1/48`); `W4` (7.1% gauge-projected transverse negative, SU2 crossover); `PROOF-EDGE-E2`, `G2G-TRANSFER-SPECTRAL` §4 (bottom-blind two-level `λ₁−λ₂`); `SIGNED-FOURIER-BLIND` §0 (positivity tools blind to atoms); `TRANSPORT-BLIND` (§3a spectral-gap≠atom, identical failure); `docs/PROOF-SKELETON.md` §0,§3 (single species = variance upper bound); `docs/YANG-MILLS-MASTER.md` lines 117,145. External (pre-vetted in-repo): Cheeger, in *Problems in Analysis* (1970) 195; Buser, *Ann. Sci. ENS* **15** (1982) 213; Bérard–Besson–Gallot, *Ann. Scuola Norm. Sup. Pisa* **9** (1982) 341; Gromov–Lévy; Ledoux; Lott–Villani, *Ann. of Math.* **169** (2009) 903; Sturm, *Acta Math.* **196** (2006) 65,133; Savaré, *Discrete Contin. Dyn. Syst.* **34** (2014) 1641; Sinclair, *Combin. Probab. Comput.* **1** (1992) 351; Jerrum–Sinclair, *SIAM J. Comput.* **18** (1989) 1149; Doyle–Snell, Carus **22** (1984); Feller, *Ann. of Math.* **55** (1952) 468; Wiener, *J. Math. Phys.* **3** (1924) 127; Shen–Zhu–Zhu, *CMP* **400** (2023) 805, arXiv:2204.12737; Savvidy, *PLB* **71** (1977) 133.
