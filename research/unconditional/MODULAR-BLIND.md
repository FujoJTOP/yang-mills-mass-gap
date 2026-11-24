# MODULAR-BLIND — independent blind red-team on the "Tomita–Takesaki modular theory / operator-algebraic QFT" species for the 4D pure SU(N) mass gap

> 2026-11-23 · **判词：SPECIES DEAD — and ALREADY-BANKED (triple-banked). This is not a 10th species; it is the modular/AQFT face already blinded and banked as `MODULAR-TAKESAKI`, `AQFT-CAUSAL`, and `H1-T2-AQFT`.** The proposal claims modular structure (modular Hamiltonian `K_O=−log Δ_O`, Bisognano–Wichmann (BW) geometric flow, **deviation from geometric = mass**) is *algebra/KMS, a different category from measure positivity*, and so might lower-bound the mass where positivity failed (anti-localization, `TRANSPORT-CORE` §4). The "different category" is real but is still **positivity** (Δ>0, KMS, strong subadditivity, monotonicity of relative entropy are all positivity/entropy statements), and it is **anti-localized in the *worst* direction** — modular data is UV/area-law (the entangling surface), the gap is IR (`k→0`). **All four kill targets land, three of them independently fatal.** Kill 1 is sharpened to airtight: local algebras are Type III₁ (Fredenhagen 1985), so `spectrum(Δ_O)=[0,∞)` and `K_O` has **no spectral gap for massive *or* massless** — the "modular spectral gap" is identically absent, so "mass = modular gap" is false (it reads `0` for every `m`). Retreating to "deviation of `K_O` from the geometric boost = mass" is circular (Kill 2): the modular data is a functional of the state `ω_m`, whose mass content is the two-point function, whose IR decay rate *is* `m` = the gap. Modulo algebra language the only mass-sensitive datum is `inf supp ν=m` = H1 (Kill 4). **Survival score 1/5** (only nonperturbativity; identical to `AQFT-CAUSAL`). Exact stuck step: a same-channel, delocalized, **mass-sensitive** modular variance bound `Var_ω(a_{k→0})≤1/λ` — which does not exist; every proven modular inequality is mass-blind, wrong-sign, or wrong-channel.

---

## §0. The tool, stated precisely (kill the real thing, not a strawman)

A Haag–Kastler net `O↦𝒜(O)` (von Neumann algebras of observables localized in region `O`) with vacuum vector `Ω`. Tomita–Takesaki supplies, for each `O` with `Ω` cyclic/separating on `𝒜(O)`: the modular operator `Δ_O`, modular conjugation `J_O`, modular group `σ_t=Ad Δ_O^{it}`, modular Hamiltonian `K_O=−log Δ_O`, and the KMS condition at `β=2π`. Bisognano–Wichmann: for a wedge `W`, `Δ_W^{it}=U(Λ_W(−2πt))` (the Lorentz boost), `J_W=CPT·R`. The proposal's target is the same object every species targets: a sign-controlling lower bound

```
target:  m_eff²(k→0) ≥ λ>0   ⟺   Var_μ(a_{k→0}) ≤ 1/λ   ⟺   atom f₁²>0   (H1 / H*; PROOF-SKELETON §0)
```

with the **novel** claims: **(A)** modular theory is a *different category* from measure positivity (algebra/KMS vs RP/LSI/FKG/Prékopa/Bakry–Émery), so it may act where positivity is blind; **(B)** the *deviation of the modular flow from the geometric (boost/conformal) flow* is a mass-sensitive functional, and (appealing to `TRANSPORT-CORE` §4 anti-localization) it might be largest exactly where the gap lives. §1–§4 are the four kill targets; §4.5 attacks claim (A)+(B) head-on. Each maps to the five constraints: (1) nonperturbative; (2) sign-controlling **lower** bound `m_eff²(k→0)≥λ`; (3) non-circular; (4) fails the free field (`f₁²=0`,`G″=0`); (5) `(L,a)`-uniform, reaches `β≈2.5`.

**Step 0 (prior version).** The exact tool is already in-repo, blinded and banked three times: `MODULAR-TAKESAKI` (Tomita–Takesaki, BW, relative entropy; two walls: mass-blind BW boost, interacting ES = gap), `AQFT-CAUSAL` (wedge modular flow / RS / split; "**refuted, not merely unproven**," by the free massless net), `H1-T2-AQFT` (operator-algebraic AQFT = "sharpest statement layer, not an independent lever"). The attack below is on *that* prior version plus this proposal's two added wrinkles (A)/(B), with the corpus's own numbers.

---

## §1. KILL 1 — Type III₁ gauges out: there is **never** a modular spectral gap (`spectrum Δ_O = [0,∞)` for **both** `m=0` and `m>0`)

**Lands — decisive, and airtight for all regions (not just wedges).**

The literal reading of "modular gap = mass gap" is a gap in the spectrum of the modular Hamiltonian `K_O=−log Δ_O`. That reading is dead before the physics starts:

- **Type III₁ is a theorem.** The local algebras `𝒜(O)` of a QFT satisfying the Wightman axioms are (isomorphic to) the hyperfinite **Type III₁** factor (Fredenhagen 1985; Yngvason review). Type is a **UV / short-distance** property of the net + state — it is *independent of the mass gap* (an IR property).
- **Type III₁ ⇒ full modular spectrum.** For a III₁ factor, the modular spectrum `S(𝒜)=[0,∞)`; for any faithful normal state (the vacuum) `spectrum(Δ_Ω)⊇S(𝒜)` and `spectrum(Δ_Ω)⊆[0,∞)` always, hence

```
spectrum(Δ_O) = [0, ∞)   ⇒   spectrum(K_O) = spectrum(−log Δ_O) = ℝ  (no gap),   ∀ O,  ∀ m.
```

So the modular Hamiltonian `K_O` is **spectrum-continuous for massive and massless alike**: its "modular gap" is identically `0`, returning no `λ>0`. (This sharpens `MODULAR-TAKESAKI` §3, which proved it only on wedges via the BW boost `K_W=U(Λ_W)`; the Type III₁ statement covers **every** bounded region.) "Mass as a modular gap" is therefore simply false — the object is mass-blind by construction.

**Verdict on Kill 1: LANDS.** Beat #2 and #5 FAIL on the spectral reading. Honest caveat: for 4D YM the net is not constructed (`G1` wall; Kill 3), so III₁ is *conditional on existence* — but conditional on existence it is the rigorous generic answer *and*, being UV, carries **zero** IR/mass information, so it cannot rescue the species.

---

## §2. KILL 2 — circular via the "deviation from geometric = mass" retreat

**Lands — this is the horn that survives even after the defender abandons "spectral gap."**

The defender retreats: "fine, the *spectrum* is continuous — but `K_O − K_O^{geom}` (deviation of the modular Hamiltonian of a bounded region from the geometric boost/conformal generator) is mass-sensitive; bind *that* below." Two independent closures:

**(a) It is a functional of the state, and the state's mass content *is* the gap.** The modular data `(Δ_ω, J_ω, K_ω)` is determined by the state `ω` and the algebra. The state's mass content is the two-point function `W_ω(x−y)=⟨Ω|φ(x)φ(y)|Ω⟩`, whose IR decay rate is the mass: `m = −lim_{|x|→∞} ∂_{|x|} log|W| = inf supp ν`. So "extract a *lower* bound on `m` from modular data" is "extract the IR decay rate of the two-point function" — which is the clustering rate — which is `m` — which is H1. Reading `m` off `Δ_ω` is reading `m` off `ω`, and `ω`'s only mass-sensitive IR datum is the gap itself.

**(b) The free-field deviation is literally `∝ m²` — quantitative only *given* `m`.** For a free massive scalar, the finite-region modular Hamiltonian deviates from the CFT one with leading correction `δK_O ∝ m²∫_O (…)+O(m⁴)` (the standard small-mass / small-region modular-Hamiltonian expansion; Casini–Huerta). Its size is a *function of* `m`: to prove `‖δK_O‖≥δ` you must already know `m≥√(δ/coeff)`, and the coefficient is theory-dependent. The "deviation" *measures* `m`; it does not lower-bound it from independent data. And for **interacting 4D YM** no such expansion exists at all (Kill 3) — so the deviation is quantitative precisely where it is circular (free field), and unavailable precisely where it is needed (YM).

**Verdict on Kill 2: LANDS.** Beat #3 FAIL. The "deviation = mass" is the gap read off the state's two-point function, with `m` as the (unknown) input; converting `deviation≥δ` into `m≥λ` presupposes `m~λ` = the answer.

---

## §3. KILL 3 — no local-algebra control: the BW/interacting-net input is itself unproven (bonus kill)

**Lands — but it is *not* the load-bearing kill; report the ordering honestly.**

- **No 4D YM net exists.** Modular theory has nothing to act on until the Haag–Kastler net of 4D pure YM is constructed — the `G1` existence wall (PROOF-SKELETON §2; `H1-T2-AQFT` §3).
- **BW for *interacting* massive 4D theories is not proven.** Bisognano–Wichmann is a theorem for **free** fields (BW 1975/76) and for conformal fields; the wedge-modular-flow = boost statement for general interacting massive Wightman fields is a conjecture (well-supported via the Unruh/KMS heuristic, but not proven). So the "geometric flow" reference point — against which "deviation = mass" is measured — is an **input**, unproven for the very (interacting, non-conformal) regime the gap lives in.

**Verdict on Kill 3: LANDS as a bonus.** Beat #1's *input* is conditional. But note: Kills 1/2/4 close the species **even granting** full modular theory and proven BW — Kill 3 only adds that even the grant is unavailable for 4D YM.

---

## §4. KILL 4 — reduces to spectral / H1: a faithful dialect, zero slack

**Lands — structurally identical to `TRANSPORT-BLIND` (Bakry–Émery = H1) and `CHEEGER-BLIND` (Cheeger = H1).**

Modulo the algebra language, the **only** mass-sensitive datum in the entire modular apparatus is `inf supp ν = m = ` the clustering rate. Everything else (BW boost, KMS `2π`, `J=CPT·R`, Type III₁, RS denseness, split/nuclearity) is mass-blind: true of the free **massless** net, which is a **counterexample**, not a bystander (`AQFT-CAUSAL` §2–§4). Concretely:

- The "modular Hamiltonian gap," once one strips the III₁ continuum (Kill 1), can only mean the **transfer-matrix gap** of the physical Hamiltonian — and on the coarse zero mode that is `G″(k→0)=1/χ=1/Var(a_{k→0})` = H1 verbatim (`OTHER-H1-POINCARE`; Ward `G″·Var(a)=1`).
- So "prove the modular flow has the right spectral/energy content" and "prove `Var_μ(a_{k→0})≤1/λ`" are the same sentence. It is a **rename** of H1, lossless (matches `MODULAR-TAKESAKI` §5: "a RENAME of the (CP)/H1 wall in modular language").

**Verdict on Kill 4: LANDS.** Beat #3 and #5 FAIL. Every proven modular inequality either gives an **upper** bound on local level density (split/nuclearity: GAP-G2-edge direction) or lower-bounds a **local** `⟨T_{kk}⟩` in a **gapless** CFT (QNEC: wrong sign, wrong channel, `AQFT-CAUSAL` §2, `H1-T2-AQFT` §4). None gives a same-channel, delocalized, mass-sensitive `Var_ω(a_{k→0})≤1/λ`.

---

## §4.5. The "different category / anti-localization" novelty (A)+(B) — killed on its own terms

The proposal's one genuine edge is **(A) "different category"** + **(B) "deviation is anti-localized to where the gap lives (`TRANSPORT-CORE` §4)."** Both are real observations; both are non-levers, and (B) is *backwards*:

**(i) The "different category" is still positivity.** Modular theory rests on: `Δ>0` (positivity of the modular operator), KMS positivity (the state is positive), **strong subadditivity** (Lieb–Ruskai), and **monotonicity of relative entropy** (Uhlmann). These are positivity/entropy statements in operator-algebraic dress. The corpus's load-bearing fact (`PROOF-SKELETON` §0; MASTER lines 117,145) is that **positivity gives the wrong direction**: `G″ = E[H] − Var + W`, and every positivity/entropy tool supplies a *lower* edge for the **subtracted** term `Var` — i.e. an *upper* bound on the variance (or a decay of a variance), never the needed *upper* bound on `Var` that *is* H1. Modular theory does not change this; it inherits it. The category difference (operator algebra vs measure) is **cosmetic, not structural** — both are positivity.

**(ii) It is anti-localized in the *worst* direction — worse than Bakry–Émery.** `TRANSPORT-CORE` §4's wall is: every positive compensation in the `Γ₂`/measure decomposition is largest at high `k`, high field, or the log-cut — i.e. **away** from the gap cell (typical field × `k→0`). For modular theory the mass-sensitive data is even more decoupled: the modular Hamiltonian `K_O` is a **short-distance, area-law** object, formally UV-divergent and dominated by entanglement across the entangling surface `∂O` at distances `~` the cutoff. The mass gap is an **IR** object (`k→0`). So modular mass-sensitivity sits at the **UV / entangling surface**, maximally displaced from the soft IR mode where `ρ` is set — a **more severe** anti-localization than BE's high-`k` vs low-`k`. The route does not dodge `TRANSPORT-CORE` §4; it reproduces it in the opposite, worse direction. (Area-law UV domination of `S` and of `K_O` is the same reason bare entanglement entropy cannot see the `k→0` threshold: `MODULAR-TAKESAKI` §4a.)

**Verdict on (A)+(B): LANDS.** The "different category" is genuine but remains positivity (wrong direction for the variance); the "anti-localization" hope is inverted — modular data is UV/area-law, the gap is IR, so the deviation lives *furthest* from the gap cell.

---

## §5. Five-point verdict (honest)

| # | Constraint | Verdict | Exact failure |
|---|---|---|---|
| 1 | Nonperturbative | **PASS** | Tomita–Takesaki / BW / KMS / III₁ are exact; no loop, all `β`. |
| 2 | Sign-controlling **lower** bound `m_eff²(k→0)≥λ` | **FAIL** | `spectrum(Δ_O)=[0,∞)` (III₁, Fredenhagen) ⇒ modular gap `≡0` for all `m` (§1); "deviation" gives `m` back, not a bound (§2); QNEC is wrong-sign/wrong-channel. |
| 3 | Non-circular | **FAIL** | Modular data = functional of `ω_m`; `ω_m`'s mass content = two-point decay rate `=m=` H1 (§2); "modular gap" = transfer-matrix gap = H1 (§4); glueball relative entropy presupposes the glueball (=the spectrum). |
| 4 | Non-SUSY + **fails free field** | **FAIL (refuted)** | The free **massless** net satisfies *every* modular hypothesis (III₁, BW, KMS `2π`, RS, split) yet is **gapless**: "modular structure ⇒ gap" is **refuted**, not unproven (`AQFT-CAUSAL` §2). Diagnostically returns `0` at `m=0` *only vacuously* (and equally `0` for free massive — mass-blind). |
| 5 | Uniform `≥ρ`, reach `β≈2.5` | **FAIL** | No `λ` exists; the modular spectrum is the mass-independent boost continuum; nothing uniform in `(L,a)`. |

**Survival score: 1/5** — passes only nonperturbativity (1); fails (2),(3),(5), and (4) is *refuted*. Same profile as `AQFT-CAUSAL` (1/5). (Most charitable reading — counting the `m=0⇒0` diagnostic as a "vacuous free-field pass" — reads **2/5**, matching `MODULAR-TAKESAKI`; but the *new* claim under attack, "deviation = mass / different category," is the refuted incarnation, so **1/5**.)

---

## §6. Why this is a re-tread, not a 10th species (decisive prior art — already banked)

The proposal is the **modular/AQFT face**, already blinded to death and banked three times, none in a new class:

- **`MODULAR-TAKESAKI`** — the *same tool* (Tomita–Takesaki, BW, relative entropy). Verdict: acid test passes only for the free Gaussian field; **two walls**: (analytic) wedge `K_W=boost` is **mass-blind** (`spectrum` identical for `m=0` and `m>0`); (structural) the only non-vacuous form — finite-region **entanglement-spectrum saturation ⟺ gap** (Li–Haldane) — is the gap relocated and is proven only for Gaussian systems, not interacting 4D YM. "A RENAME of the (CP)/H1 wall in modular language."
- **`AQFT-CAUSAL`** — wedge modular flow / RS / split/nuclearity. Verdict: BW boost is **mass-blind**; Reeh–Schlieder is **anti-clustering** (sharpest at the gapless free field); split/nuclearity is UV (upper bound, wrong direction); the free massless net satisfies all structural hypotheses yet is gapless ⇒ the implication is **refuted**. The sole mass-sensitive datum = clustering rate = H1.
- **`H1-T2-AQFT`** — operator-algebraic AQFT as a lever. Verdict: "**sharpest statement layer, not an independent lever**"; modular structure is kinematic and mass-blind; every gapped AQFT case has the mass as **input**; no modular statement lower-bounds `f₁²` in 4D.
- **`C8-SPECIES-NOGO`** — the *empty-referent* pattern the modular route repeats: a structural class that is **true of the gapless free field** cannot separate gapped from gapless; the "modular mass certificate" has no referent that isn't the gap.

The proposal's two additions — **(A)** "different category from measure positivity" and **(B)** "deviation anti-localized to the gap cell" — are **already analyzed** (`AQFT-CAUSAL` §2 is exactly (A), refuted; `TRANSPORT-CORE` §4 is the wall (B) claims to dodge). §4.5 shows (A) is still positivity and (B) is inverted. So the "10th species" is the 4th blinding of the **same already-banked** modular/AQFT species — `MODULAR-BLIND` is a stamp, not a new kill.

---

## §7. What, if anything, survives (the honest residue)

- **Genuine, already-banked, not new:** Tomita–Takesaki + BW (`K_W=U(Λ_W)`, KMS `β=2π`, `J=CPT·R`), Type III₁ of local algebras (Fredenhagen), and the **statement layer** — AQFT's sharpest *definition* of "gap ⇔ atom in the KL measure ⇔ a mass shell exists" (Haag–Swieca; Buchholz–Fredenhagen; `H1-T2-AQFT` §4). These are real; none is a lever.
- **A faithful dialect, not a closure:** modular language *is* the cleanest way to say "the gap is a mass shell / an atom / a clustering rate" (= H1/H*). It lets one cite Haag–Ruelle, Haag–Swieca, BGL machinery. It is **lossless** — it cannot weaken H1 ("H1 ⇒ gap, H1 unreduced," `REDTEAM-FORGE` cross-cutting).
- **The exact missing OBJECT** is one **mass-sensitive** modular inequality in the `k→0` channel: `Var_ω(a_{k→0}) ≤ 1/λ` (equivalently, a same-channel coercion `T ≤ C·f₁²`) — every proven modular inequality is mass-blind (BW/KMS), upper-only (split/nuclearity), or wrong-sign/wrong-channel (QNEC). This is the same "variance upper bound" species the skeleton names at G2.5/H1 and G1.2/L3 (`PROOF-SKELETON` §0,§3; `YANG-MILLS-MASTER` lines 117,145). Modular theory does not supply it; it re-expresses the demand for it in operator-algebraic costume.

---

## §8. Sharpest line

**The modular species dies on the nature of the modular operator itself: local algebras of 4D QFT are Type III₁, so `spectrum(Δ_O)=[0,∞)` and `K_O=−log Δ_O` has a *continuous* spectrum — "no invariant vectors, no modular gap" — for massive and massless theories alike (Fredenhagen 1985; `MODULAR-TAKESAKI` §3 sharpened from wedges to all regions). The "modular gap = mass" is therefore not merely unproven but *false*: the object is mass-blind and reads `0` for every `m`. The retreat to "deviation of `K_O` from the geometric boost = mass" is circular: the modular data is a functional of `ω_m`, whose only mass-sensitive IR datum is the two-point function's decay rate `=inf supp ν=m=` the gap (`MODULAR-TAKESAKI` §3; `AQFT-CAUSAL` §5). The "different category" is still positivity — `Δ>0`, KMS, strong subadditivity, monotonicity of relative entropy — which gives the *wrong direction* for the subtracted variance `Var` (PROOF-SKELETON §0). And its "anti-localization" hope is inverted: modular data is UV/area-law (the entangling surface), the gap is IR (`k→0`), so the mass-sensitive deviation is largest *furthest* from the gap cell — a worse `TRANSPORT-CORE` §4, not an escape from it. The free massless net satisfies every modular hypothesis yet is gapless: the implication "modular structure ⇒ gap" is **refuted**, not merely open (`AQFT-CAUSAL` §2). Modular theory is the sharpest *language*, the weakest *lever*. Correct language, same wall — and a wall already banked.**

**Exact stuck step.** A same-channel, delocalized, **mass-sensitive** modular variance bound
```
∃λ>0:  Var_ω(a_{k→0}) ≤ 1/λ   [edge atom, not local energy density]   — uniform (L,a),
```
equivalently a coercion `T ≤ C·f₁²` (= H* at `r₀>0`). Every known modular inequality is mass-blind (BW/KMS), upper-only on local level density (split/nuclearity), or lower-bounds the wrong channel in a gapless theory (QNEC). Producing it *is* H1/H*; modular theory does not produce it, it *is* the demand for it.

---

*Refs.* Internal (read/cited as ground, not invented): `MODULAR-TAKESAKI` (Tomita–Takesaki/BW/relative entropy; wedge boost mass-blind; interacting ES = gap = H1); `AQFT-CAUSAL` (wedge modular flow mass-blind; RS anti-clustering; split/nuclearity UV-upper; free massless net refutes "structure⇒gap"; mass-sensitive datum = clustering rate = H1; QNEC wrong-sign/wrong-channel); `H1-T2-AQFT` (AQFT = statement layer, not lever; every gapped AQFT case has mass as input); `C8-SPECIES-NOGO` (empty-referent: mass-blind class cannot separate gapped/gapless); `TRANSPORT-CORE` §4 (anti-localization of every positive compensation); `TRANSPORT-BLIND` (BE = H1 faithful dialect); `CHEEGER-BLIND` (Cheeger = H1); PROOF-JACOBIAN-SAVES / ABC-NEW-INEQ (`G″=E[H]−Var+W`; Ward `G″·Var(a)=1`); OTHER-H1-POINCARE (coarse `G″=1/χ`); GAP-G2-edge (split/nuclearity ⇒ upper only); `docs/PROOF-SKELETON.md` §0,§3 (single species = variance upper bound, UV=L3/IR=H1); `docs/YANG-MILLS-MASTER.md` lines 117,145; G1-ASSEMBLE-100 (no 4D YM net). External (pre-vetted in-repo): Haag–Kastler, JMP **5** (1964) 848; Bisognano–Wichmann, JMP **16** (1975) 985 / **17** (1976) 303; Tomita (1967) / Takesaki (1970); **Fredenhagen, CMP 97 (1985) 79** (local algebras Type III₁); Yngvason, *The role of type III factors in QFT*, Rept. Math. Phys. **55** (2005) 135; Araki, PTP **32** (1964) 844; Haag–Hugenholtz–Winnink, CMP **5** (1967) 215; Reeh–Schlieder, Nuovo Cim. **22** (1961) 1051; Borchers, CMP **143** (1992) 327; Wiesbrock, CMP **157** (1993) 83; Brunetti–Guido–Longo, RMP **14** (2002) 759; Buchholz–Wichmann, CMP **106** (1986) 321; Summers–Werner, CMP **110** (1987) 247; Haag–Swieca, CMP **1** (1965) 308; Buchholz–Fredenhagen, CMP **84** (1982) 1; Li–Haldane, PRL **101** (2008) 010504; Casini–Huerta, J. Phys. A **42** (2009) 504007 (finite-region modular Hamiltonian; small-mass correction `∝m²`); Balakrishnan–Faulkner–Khandker–Wang, JHEP **09** (2019) 020 (QNEC); Ceyhan–Faulkner, arXiv:1812.04683; Lieb–Ruskai (strong subadditivity); Uhlmann (monotonicity of relative entropy).
