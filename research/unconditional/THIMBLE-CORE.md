# THIMBLE-CORE — candidate 7th species, Picard–Lefschetz thimble / resurgence across the Savvidy crossover: the phase is *not* annihilated by `f₁²=|A₁|²` (it survives as inter-thimble interference), and that is precisely why the species dies — what survives is the undetermined Stokes constant = the gluon condensate = the atom

> 2026-11-23 · Core adjudication of the seventh species: *a rigorous thimble decomposition of the YM path integral across the Savvidy crossover, whose genuine PHASE `e^{i Im S_σ}` (not a positivity) lower-bounds the atom `f₁²=|⟨0|trF²|G⟩|²`.* **判词 (verdict): DEAD as a lower-bound lever; score 1.5/5 — the highest-quality DEAD on the board.** It is the only species that is genuinely *non*-positivity (oscillatory/complex-analytic, not a real Schur-concave inequality) and it is the only one that is *free-field-safe by silence rather than by vacuity* (§3, condition 4 PASS). But the crux resolves the wrong way for a bound, on two distinct rails. **(Rail I — global phase.)** If the lever is a single thimble's phase, `f₁²=|A₁|²` annihilates it — verbatim the death of species 6. **(Rail II — relative/interference phase.)** This one *does survive* the modulus-square: distinct thimble contributions to the *same* amplitude `A₁` interfere, `|Σ_σ A_σ|² = Σ|A_σ|² + 2Re Σ_{σ<σ'} A_σ \bar A_{σ'}`, and `Re(A_σ\bar A_{σ'})` carries `cos(Im S_σ−Im S_{σ'})`. So the thimble species is the only candidate yet whose "phase" is *not* deleted by `|·|²`. **But the interference term is (a) sign-indefinite — it can cancel the diagonal as easily as reinforce it, i.e. one-sidedness again — and (b) governed by the relative phase/Stokes datum `σ_S`, which on ℝ⁴ is a free trans-series parameter (no bion/monopole sector; Maiezza–Vasquez 2405.01639), equal to the gluon condensate `⟨TrF²⟩`.** So the abundance that survives the modulus-square is exactly the abundance that is undetermined, and the undetermined datum is the atom. The thimble sum is a formal trans-series with free Stokes constants, *not* a rigorous lower bound; it reduces to the same wall. Exact stuck step unchanged: `σ_S` located and sign-definite = `⟨TrF²⟩` = `f₁²` = H1.

---

## §1 · The thimble decomposition, written out — and exactly where a lower bound could live

**1.1 Complexified config space and a holomorphic action.** Take links `U_e∈SL(N,ℂ)` (or `A_μ∈sl(N,ℂ)` in the continuum) and the analytically continued action
`S_θ(A) = (1/2g²)∫ Tr F_μν F_μν + iθ(1/8π²)∫ Tr F F̃`, `F=dA+A∧A`.
`S_θ` is holomorphic in `A` (no `\bar A`), so the Morse–Picard–Lefschetz machine applies to `Z=∫_{Γ} e^{−S_θ}𝒟A` after the (unproven for 4D gauge) deformation of `Γ` from the real contour (Witten 1001.2933).

**1.2 Saddles `∂S=0 ⇔ D_μF^{μν}=0` (complexified).** Three physical families:
- **(i) Perturbative** `A=0`, `S=0`, `Hess=−D²` PSD (kernel = constant gauge orbits) — the Gaussian/one-loop vacuum.
- **(ii) BPST (anti)instantons**, (anti)self-dual ⇒ BPS ⇒ *global minima*, `S=8π²|Q|/g²`, `Hess` PSD on the quotient up to moduli (zero modes). Stable thimbles; phase trivial (`Im S=0`).
- **(iii) Savvidy / Nielsen–Olesen family** — constant covariantly-constant abelian `F=B` in a Cartan direction (a *cone*: the Cartan direction and `B` are moduli, `B` fixed at one loop by `V_eff`). `S_cl=(B²/2g²)V`, volume-extensive and real. The Hessian on transverse fluctuations is the Landau ladder
  `H=−D²+2\,ad\,F`, `λ_{n,s}=k_∥²+(2n+1)g|B|+2s\,gB`, bottom `λ_min^{}=-g|B|` (PROOF-SAVVIDY-EXACT-SIGN). `N_neg/V∝(gB)²`; on the lattice the gauge-projected transverse Hessian is genuinely indefinite (7.1% negative, min `−0.82β` at the SU(2) crossover, W4-blocked-positivity; full microscopic min `−1.12β`, 22–49% negative at SU(3), W6-BE-at-2.5). **These are the saddles with `N_neg=1` that generate the thimble phase.**

**1.3 The decomposition and its measure.** `Z=Σ_σ n_σ∫_{J_σ}e^{−S_θ}𝒟A`, `J_σ` = downward-flow (steepest-descent) thimble, `n_σ∈ℤ` = intersection numbers of the original contour with `J_σ`. On `J_σ`, `Im S_θ≡Im S_θ(φ_σ)` (the flow conserves it), so
`∫_{J_σ}e^{−S}dμ = e^{−i Im S_σ(φ_σ)} e^{−i θ Q_σ} ∫_{J_σ} e^{−Re S}\,ω_σ\,dμ_σ`, `ω_σ=\arg\det J_σ` (the residual Jacobian phase; PROOF-THIMBLE).

**1.4 One-loop determinants and the phase.** Gaussian around each saddle:
`A_σ = e^{−S_σ}\,\big(\det{}' H_σ\big)^{−1/2}\,(i)^{N_neg}\,(\text{collective/zero-mode});`
each negative eigenvalue contributes a factor `±i/√|λ_-|` (the steepest-descent contour over the unstable direction is rotated to the imaginary axis); the collective zero modes (gauge orbits — **Morse–Bott**, non-isolated; Gribov copies on any slice, Gribov 1978) give the Faddeev–Popov/measure factors with a sign that flips across the Gribov horizon. For the Savvidy saddle `N_neg=1` ⇒ the leading phase is the one-loop Nielsen–Olesen `i`-factor, `Im V_eff∝(gB)²`.

**1.5 Stokes phenomena.** As `β,θ,L` vary, saddles enter/leave the set with nonzero `n_σ` across **Stokes lines**; the thimble sum jumps by `σ_S·A_{σ'}` with **Stokes constant `σ_S`** (Berry, Proc. R. Soc. A **422** (1989) 7; Écalle; Costin 2009; Aniceto–Başar–Schiappa, Phys. Rept. **809** (2019)). These `σ_S` are the ambiguities; on ℝ⁴ no monopole/bion sector fixes them (Maiezza–Vasquez, PLB 2024, 2405.01639).

**1.6 Exactly where a lower bound could live — the three candidate sites.** The atom `f₁²=|A₁|²`, `A₁=Σ_σ A_σ` (thimble contributions to the operator amplitude `⟨0|trF²|G⟩`). Expanding:
```
f₁² = Σ_σ |A_σ|²  +  2 Re Σ_{σ<σ'} A_σ \bar A_{σ'}          (diagonal + interference)
```
- **(S-i) Diagonal (magnitude) floor.** `|A_σ|² = e^{−2 Re S_σ}\,|det{}' H_σ|^{−1}(⋯)²`. The exponential is *real positive* — the only positivity-shaped piece in the whole species — and it *vanishes at `g=0`* (good, free-field). A floor `≥ρ_σ>0` would need `Re S_σ` bounded above and `|det{}'H_σ|` bounded — both are the *dynamical* background content, and the overall normalization is the trans-series/Stokes constant (`§2`). **Reduces to `σ_S`.**
- **(S-ii) Interference.** `2Re(A_σ\bar A_{σ'})` carries `cos(Im S_σ−Im S_{σ'})`; it is **sign-indefinite** (destructive cancellation possible) and the phase is the undetermined Stokes datum. A *lower* bound needs constructive non-cancellation — *one-sidedness* in interference clothing. **Reduces to `σ_S` + the wall.**
- **(S-iii) Resummed (median / lateral Borel) sum.** The resummed sum is `g`-continuous and `0↦0` at `g=0` (CONTINUUM-ASYMPTOTIC-BEYOND §4); a homomorphism of the ring fixed on `0` cannot manufacture `f₁²≥r₀>0`. **Empty by shape.**
So the decomposition has exactly one positivity-shaped piece (S-i), and it is already the condensate at `g=O(1)`.

---

## §2 · The phase-vs-modulus-square question — resolved on two rails, both dead for a bound

This is the crux where species 6 died. The thimble species forces it to *two* sub-questions, and it answers them differently — which is the whole reason it is the correct 7th species to test and the whole reason it still dies.

**Rail I — a single (global) thimble phase is annihilated.** If the lever is "the dominant thimble carries a relative sign / global phase `e^{iφ}`," then `A₁→e^{iφ}A₁` and `f₁²=|A₁|²` is *invariant*. This is the literal death of species 6 (SIGNED-FOURIER-CORE `§1`; SIGNED-FOURIER-BLIND §1a: "a modulus-square annihilates every phase"). **Rail I dies identically.** Any appeal to "the Savvidy thimble is `i` relative to the perturbative one" as a *sign* for `f₁²` is this rail and is dead.

**Rail II — distinct thimbles interfere, and the interference survives `|·|²`.** This is the one place the thimble species differs *structurally* from every other species tried. If there are ≥2 saddle contributions to the *same* amplitude `A₁` (perturbative `A_0` + Savvidy `A_\text{sav}` + …), then
```
f₁² = |A_0 + A_\text{sav} e^{iφ} + …|² = |A_0|² + |A_\text{sav}|² + 2|A_0||A_\text{sav}|\cos φ + …
```
and `cos φ`, `φ=Im S_\text{sav}−Im S_0`, is **not** deleted by the modulus-square. The phase is real, visible, and (unlike species 6) cannot be removed by `A₁→−A₁`. **So the honest answer to "does the thimble phase survive the modulus-square?" is: as a global phase, no (Rail I); as a relative/interference phase, yes (Rail II).**

**But Rail II yields no lower bound, for two independent reasons.**
1. **Undetermined.** `φ` (equivalently the relative Stokes constant `σ_S` between the perturbative and Savvidy thimbles) is free on ℝ⁴ (no bion; RESURGENCE `§2`; Maiezza–Vasquez 2405.01639). The interference term is a real but *unfixed* number.
2. **Sign-indefinite (one-sidedness).** `2|A_0||A_\text{sav}|\cos φ` is bounded *below by* `−2|A_0||A_\text{sav}|`; a floor on `f₁²` requires the interference not to over-cancel — i.e. `cos φ ≥ −1+δ` — which is precisely a **coherence / one-sidedness** statement, the same wall that killed the six positivity species (SCHUR-CONVEX: positivity is Schur-concave and cannot force an atom; here the interference is an *indefinite* correction, the phase-space version of that obstruction).
3. **Unproven to exist.** Even to *write* `n_\text{sav}≠0` you need the thimble decomposition equality at `β≈2.5`, unproven outside `β<β₁` (T5; PROOF-THIMBLE).

**§2 outcome.** The modulus-square does *not* trivially annihilate the thimble phase — it survives as interference. That is a genuine, new structural fact, and it is why thimbles were the right 7th species to test. But the surviving object is the interference term, whose control is exactly the Stokes datum = the gluon condensate = the atom. Rail I reduces to species 6; Rail II reduces to the Borel ambiguity. **No rail gives a rigorous lower bound; both give "the series is fixed only up to exactly the physics it omits."**

---

## §3 · The two acid tests — outcomes

### Acid test 1 — without a measure-positivity inequality **and** non-circular in the condensate/OPE? **PARTIAL (passes "not positivity", FAILS "non-circular") ⇒ DEAD.**
- **Not positivity: PASS (the genuine virtue).** Unlike all six prior species, the thimble sum is not a Schur-concave variance/decay bound; it is an oscillatory complex-analytic decomposition. It is *structurally* capable of carrying the atom (unlike RP/LSI/FKG/BL/Bochner/entropy). This is the one real reason it is the seventh species.
- **Non-circular: FAIL.** The magnitude that `f₁²` would extract is fixed by the Stokes constant `σ_S`, and by the renormalon factorisation `σ_S = ⟨TrF²⟩` (dim-4 `n=1`; CONTINUUM-ASYMPTOTIC-BEYOND §1–2, RESURGENCE §2). To bound `f₁²` below you must *fix* `σ_S`, i.e. the gluon condensate = the atom. **Circular.**
- **The deeper miss — there is no bare sign problem to escape.** At real `β≈2.5` the physical weight `e^{−S_W}>0` and is reflection-positive (PROOF-THIMBLE premise check); PL is machinery for *complex/oscillatory* weights. So the "escape from positivity" is illusory: by **Cauchy invariance** any valid thimble contour reproduces the *same* positive-RP correlators, hence the *same* gap theorem at `β<β₀` and *nothing* at `β≈2.5`. A contour deformation cannot add convergence the physical theory lacks (PROOF-THIMBLE; W4-cert-vs-gap: the crossover is one *sufficient* route, necessity is false).
**Outcome: DEAD.** The species is not positivity (true), but it is circular in the condensate **and** it deforms an integral that has no sign problem, so it returns the same theory.

### Acid test 2 — does it VANISH on the free field (`g→0`, no Savvidy instability, trivial thimble, `f₁²=0`)? **PASS — but by silence, not by the thimble.**
- As `g→0`: the classical action of any nonzero saddle is `(1/2g²)∫TrF²→∞`, so only the perturbative saddle `A=0` survives; the Savvidy cone collapses (`B_min∝Λ e^{−c/g²}→0`) and the negative-mode `(i)`-factor disappears (`N_neg/V∝(gB)²→0`). The thimble sum returns the Gaussian one-loop vacuum.
- There `A₁=⟨0|trF²|\text{single gluon}⟩=0`: `trF²=F\!\cdot\!F` is a **two-gluon** operator with no one-particle overlap (the free KL measure is a two-photon cut, `f₁²=0`; PROOF-F1SQ). So `f₁²→0` correctly. **It does not return `>0` on the free field ⇒ not vacuous** (contrast species 6, whose benign Bochner face held *strictly* on the free heat kernel). **Acid test 2 PASSES.**
- **Honesty caveat.** The pass is *silent*: the vanishing is carried by the operator identification (`A₁=0` because `trF²` is 2-particle), which is true of *every* species, not by anything specific to thimbles. The thimble structure is *neutral* on the free field — it just computes the Gaussian. So the discriminator that separates interacting (`f₁²>0`) from free (`f₁²=0`) is imported (= the atom), not produced. This is the same "discrimination imported, not earned" issue that made species 6's condition 4 *hollow* on the strict read — except here the tool returns the *correct* `0` instead of a *wrong* `>0`, which is why it earns the box (§5 condition 4) rather than failing it.

---

## §4 · The crossover attempt (`β≈2.5`) — and the three forge routes (charter rule 1)

At the crossover the Savvidy saddles proliferate across Stokes lines; the negative mode `λ_min=−g|B|` makes one-loop `V_eff` complex, `Im V_eff∝(gB)²` — the leading **Borel singularity = IR renormalon** at `u=2n/β₀` on the positive axis, factorising into `⟨TrF²⟩` `(n=1)`. The two lateral Borel sums differ by the trans-series amount `2πiΣ_n C_n e^{−2n/(β₀g²)}`, `C_n∝γ_n⟨O_{4n}⟩` (CONTINUUM-ASYMPTOTIC-BEYOND §1–2). The thimble sum across the crossover *is* this resurgent trans-series.

**Attempt to read the atom off it.** The glueball pole `m₁` and residue `f₁²` would be features of the Borel-resummed correlator. But the resummation is ambiguous by `σ_S`, and the ambiguity is **largest exactly at the IR edge / `k→0`**, where the `(CP)`/H1 bound bites (CONTINUUM-ASYMPTOTIC-BEYOND §4; RESURGENCE §4). The trans-series *localises* the wall to `σ_S` and sets the exponential scale but cannot fix its sign or its value. Three forge routes, each landing on the wall:
- **(α) Interference floor.** `f₁² ≥ (|A_0|−|A_\text{sav}|)²₊`. Needs `|A_0|>0` (the atom) *and* `|A_\text{sav}|` controlled; both are `σ_S`-ambiguous ⇒ reduces to H1 (the `(S-i)+(S-ii)` of §1.6).
- **(β) Analytic structure / isolation `I2`.** The thimble sum's analytic structure *is* the trans-series; its Borel singularities are the renormalons, not (generically) the glueball pole. To make the pole a stable feature of the resummed sum you must fix the ambiguous data, i.e. `σ_S` = condensate = the atom. The residue-level gap `I2` is *not* produced (PROOF-P1P2 §3, PROOF-P3 §5: grading is amplitude, not eigenvalue spacing).
- **(γ) Cauchy-invariance transport (β<1/48 → 2.5).** By Cauchy invariance the thimble reproduces the identical positive-RP correlators; so a thimble "proof" would assert the gap is β-independent across `β₁`, which is false in general (phase transitions) and, where true, is the unproven single-phase `(CP)` statement (MASTER §6). Empty.
**No forge route survives; (α)/(β) reduce to `σ_S`=H1, (γ) reduces to `(CP)`.**

*> Precision note (charter §2).* RESURGENCE `§2` writes `S_sav=8π²/g²`; that is the **BPST instanton** action, not the Savvidy classical action (which is `(B²/2g²)V`, volume-extensive) nor the renormalon order (which is `e^{−2n/(β₀g²)}=(Λ/μ)^{2n}`). The physical conclusion is unaffected; the exponential *order* is the renormalon power, and the *scale* is the Savvidy/condensate scale.

---

## §5 · Scorecard (0–5), sharpest line, exact stuck step

| # | Condition | Verdict | Why |
|---|---|---|---|
| 1 | **Nonperturbative-rigorous** | **PARTIAL (½)** | Trans-series/saddles are genuinely nonperturbative (no `1/N`, no loops) — but the 4D-gauge thimble decomposition `Z=Σn_σ∫_{J_σ}` is **unproven** off `β<β₁` (T5), the physical real-β measure has **no sign problem** (PROOF-THIMBLE), and Morse–Bott/Gribov defeat a clean single-thimble statement. Rigorous *local* Borel-singularity analysis (Costin/Écalle) does not extend to a rigorous 4D gauge *bound*. |
| 2 | **Sign-controlling LOWER bound `m_eff²(k→0)≥λ`** | **FAIL** | The phase is a phase; the surviving interference term is *sign-indefinite* and `σ_S`-undetermined. No `λ`. |
| 3 | **Non-circular** | **FAIL** | `σ_S=⟨TrF²⟩=f₁²`; fixing the Stokes constant to bound the atom *is* assuming the atom. |
| 4 | **Non-SUSY + fails free field** | **PASS** | Non-SUSY: thimbles need no BPS ✓. Free field (`g→0`): trivial thimble, `A₁=⟨0|trF²|1 gluon⟩=0` ⇒ `f₁²→0` ✓ (silent, not vacuous — the one clean pass). |
| 5 | **`(L,a)`-uniform `≥λ`** | **FAIL** | Cauchy invariance ⇒ only the `β<β₀` theorem, nothing at `2.5`; no `a→0`; Stokes lines proliferate and the ambiguity is *maximal* at the `k→0`/IR edge exactly where the bound bites. |

**Score: 1.5/5** (½×#1 + 1×#4). Reconciles with RESURGENCE's 2/5 as follows: RESURGENCE scored #1 a full PASS reading the condition as "genuinely nonperturbative"; the thimble-CORE splits off the "rigorous 4D decomposition / no-sign-problem" half, which fails — hence 1.5 rather than 2. This is the **highest-quality DEAD** on the board: it is the only species that is *both* genuinely non-positivity *and* free-field-safe, i.e. the only one not killed by the elementary modulus-square/vacuity argument *alone* — it is killed instead by the Stokes ambiguity being the very object it must compute.

**Sharpest line.** *The thimble species is the only adversary that escapes the positivity wall, and that is exactly why it fails: freed from positivity, its phase is not deleted by the modulus-square (distinct thimbles interfere: `|Σ_σ A_σ|²=Σ|A_σ|²+2ReΣ A_σ\bar A_{σ'}`), but the interference term is sign-indefinite (one-sidedness again) and governed by the relative phase — the Stokes constant `σ_S` — which on ℝ⁴ is free and equals the gluon condensate `⟨TrF²⟩`. So the one object that survives `|·|²` is precisely the one object the theory refuses to fix, and that object is the atom. Rail I (a global thimble phase) is species 6 over again; Rail II (interference) is the Borel ambiguity in thimble coordinates. Either way: "the series is fixed only up to exactly the physics it omits."*

**Exact stuck step (unchanged wall, new coordinate).**
```
Need:  a Stokes constant σ_S (the relative phase/coefficient between the perturbative A=0
       thimble and the Savvidy/Nielsen–Olesen thimble) that is
   (a) LOCATED and SIGN-DEFINITE   — not a free trans-series parameter, and
   (b) provably = the physical path-integral trans-series at β≈2.5   (thimble decomposition = measure),
   so that  f₁² = |A_pert + σ_S e^{iφ} e^{−S_sav}|²  ≥ ρ₀ > 0,  uniformly in (L,a),  → 0 at g=0.
Reality on ℝ⁴:  σ_S is FREE (no bion/monopole sector; Maiezza–Vasquez);  σ_S = ⟨TrF²⟩ (n=1 renormalon);
       the decomposition equality is unproven off β<β₁ (T5);  and the physical real-β measure has NO sign problem,
       so Cauchy invariance returns the SAME β<β₀ gap and nothing at 2.5.
⟺  σ_S is the gluon condensate  ⟺  f₁²  ⟺  H1.       Species DEAD; the wall is unchanged.
```
The missing object is **not** a thimble lemma. Thimbles are fully understood (Picard–Lefschetz; Berry/Écalle Stokes; Costin). What is missing is the **value and sign of the Stokes constant** — a dynamical/nonperturbative number identical to the gluon condensate and to the atom. Positivity fixed the six species' obstruction at the *inequality*; thimbles fix this one at the *parameter*. Same wall, and the seventh-species catalogue now records that the complex-analytic route reduces to it as cleanly as the positivity route did.

*Refs.* External: E. Witten, arXiv:1001.2933 (analytic continuation / Lefschetz thimbles); Picard–Lefschetz theory (Arnold–Gusein-Zade–Varchenko; Pham); M. V. Berry, Proc. R. Soc. Lond. A **422** (1989) 7 (Stokes phenomenon); J. Écalle, *Les fonctions résurgentes* (Orsay, 1981); O. Costin, *Asymptotics and Borel Summability* (CRC, 2009); Aniceto–Başar–Schiappa, Phys. Rept. **809** (2019) 1; Savvidy, Phys. Lett. B **71** (1977) 133; Nielsen–Olesen, Nucl. Phys. B **144** (1978) 376; Leutwyler, Nucl. Phys. B **179** (1981) 129; Gribov, Nucl. Phys. B **139** (1978) 1; 't Hooft (1977, renormalons); Beneke, Phys. Rept. **317** (1999) 1; Maiezza–Vasquez, PLB (2024) arXiv:2405.01639; Cristoforetti et al., arXiv:1205.3996, arXiv:1303.7204 (thimble sign problem). Internal: PROOF-THIMBLE (no bare sign problem; Morse–Bott/Gribov; residual Jacobian phase; Cauchy invariance — same gap); RESURGENCE (`§1` convergent below β₁ / contentless; `§2` σ free; `§3` `0↦0`); CONTINUUM-ASYMPTOTIC-BEYOND (`§1–4` renormalon = `⟨TrF²⟩`, `0↦0`, ambiguity maximal at IR edge); W4-blocked-positivity (gauge-projected transverse Hessian indefinite, 7.1% neg, min −0.82β); W4-cert-vs-gap (crossover certificate sufficient, not necessary); PROOF-SAVVIDY-EXACT-SIGN (Landau ladder `λ_min=−g|B|`, `N_neg∝(gB)²`); PROOF-SAVVIDY-WGEO-DOMINATION (negativity physical, survives quotient); OTHER-SAVVIDY-SIGN (parity ⇒ balanced Q=0 mixture = max excess); W6-BE-at-2.5 (min −1.12β, 22–49% neg); SIGNED-FOURIER-CORE / SIGNED-FOURIER-BLIND / SIGNED-FOURIER-BUILD (species-6 modulus-square kill — Rail I); C8-SPECIES-NOGO (frustration invisible; sign decoupled from bounded quantity); SCHUR-CONVEX (positivity Schur-concave; atom not forced); PROOF-P1P2 / PROOF-P3-UNCONDITIONAL (I1 ≠ I2; provable only β<1/48); PROOF-F1SQ (free `:φ²:` `f₁²=0`); MASTER §5–6 (six immunity roots; residual object = volume-uniform `m_{0⁺⁺}` lower bound = H1).
