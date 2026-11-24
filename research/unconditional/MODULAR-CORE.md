# MODULAR-CORE — 10th-species CORE attack: Tomita–Takesaki / operator-algebraic QFT. **STEP-0 finding: modular theory is already banked three times** (`MODULAR-TAKESAKI` 1/5, `H1-T2-AQFT` 1/5, `AQFT-CAUSAL` 1/5). The one genuinely-untested core — *does the **Type III₁ spectral structure** itself force a clean no-go, and does the modular route **evade anti-localization** on the sign-bearing cell?* — is decided **NO on both**: the modular spectrum `σ(Δ_W)=[0,∞)` is mass-independent by the **type** (a theory-independent statement, stronger than the banked "free boost is mass-blind"), and the modular Hamiltonian `K_W=2π∫x¹T₀₀` is a **local stress-tensor integral that vanishes on the delocalized `k→0` sign cell** — modular theory **reproduces** anti-localization rather than escaping it. **判词 (verdict): NO NEW LEVER — a clean NEW NO-GO. Score 1/5.**

> **2026-11-23 · 判词 (verdict, up front).** The vacuum Ω is cyclic+separating for each local algebra 𝒜(O); Tomita–Takesaki gives `Δ, J, K=−logΔ`, and Bisognano–Wichmann identifies the **wedge** modular flow with the Lorentz boost (`Δ_W^{it}=U(Λ_W(−2πt))`, KMS `β=2π`). This is genuine **algebra/KMS**, *not* measure-geometry positivity — so the species **passes** the first half of acid test 1 (it is not another member of the reflection-positivity/FKG family). But the deep question — *does modular structure lower-bound the mass?* — has a **structural NO** that is stronger than every banked version: local 4D algebras are (expected) **Type III₁**, whose modular operator has the **full spectrum `σ(Δ)=[0,∞)`** (equivalently `σ(K)=ℝ`), the modular flow is **ergodic** (no almost-periodic part, no invariant vector), so there is **no modular spectral gap** — and this is **mass-independent** because the type is fixed by the net+vacuum (UV/kinematic), not by `m`. Hence *any functional of the modular spectrum is a constant `0` in the gap variable*: a naive "modular-gap ⇒ mass" is not merely unproven but **structurally vacuous**. The only non-spectral modular handles are the spectral **measure** (entanglement entropy) and the **BW deviation** `K_mod−K_boost`; both are mass-sensitive **only through the correlators** — the Gaussian entanglement spectrum works (free/Gaussian only, banked `MODULAR-TAKESAKI`), and the BW deviation is an **interaction functional** whose bound flows from the low-`k` stress tensor = the KL spectral measure = `m`. The sign-bearing cell (typical field × `k→0` × `O(π)` phase, `TRANSPORT-CORE §4`) is where every structural positivity has vanished; modular theory is the **ninth costume of the same anti-localization** because `K_W` integrates the **local** energy density and the sign mode is **delocalized + soft** (`T₀₀→0` there). **Net: the honest output is a precise no-go — "Type III₁ makes the modular spectrum gapless independent of mass, and the only mass-sensitive modular data is circular" — worth stating as a theorem. 1/5.**

---

## §0. STEP-0 honesty discharge (what is already banked)

Per project rule the corpus already contains the modular / AQFT legs; this report does **not** re-derive them. It isolates the three claims none of them makes, which are the substance of a *core* attack.

- **`MODULAR-TAKESAKI.md`** (banked, 1/5): the free-field **entanglement-saturation** diagnostic. Correct acid-test pass in its only well-posed (Gaussian) form: half-chain entropy of a free-fermion chain is `S=(1/6)lnL+…` at `m=0` (log-divergent, no gap) and **saturates** (`L`-independent) at `m>0`. But the **continuum wedge** `Δ_W=boost` is mass-blind (identical for `m=0` and `m>0`), and extending "entanglement saturation ⟺ gap" to interacting 4D YM **is** Li–Haldane / the gap. Relative entropy `S(ω_glueball‖ω_vac)` is UV/area-law-dominated and presupposes the glueball (= the gap). **Establishes: the free Gaussian diagnostic works; the interacting extension is circular.**
- **`H1-T2-AQFT.md`** (banked, 1/5): AQFT is the sharpest **statement layer** (gap ⟺ KL atom ⟺ mass shell, Haag–Swieca / Buchholz–Fredenhagen), **not an independent lever**; modular structure is kinematic, mass-blind; where AQFT "proves a gap" the mass is solvable or an input (P(φ)₂, YM₂, Lechner models); no 4D YM net exists (G1). QNEC lower-bounds **local** `⟨T_kk⟩` (wrong sign *and* channel).
- **`AQFT-CAUSAL.md`** (banked, 1/5): the **modular-mass-blindness lemma** (`Δ_W^{it}=boost`, `spec Δ_W=[0,∞)`, `β=2π`, same for `m=0`/`m>0`), Reeh–Schlieder is **anti-clustering** (gapless-compatible), split/nuclearity is **UV**. Sharpest framing: the free massless net satisfies **every** structural hypothesis yet is gapless, so "modular structure ⇒ gap" is **refuted**, not merely unproven. The sole mass-sensitive datum is the clustering rate `inf supp ν=m` = H1.

**What these three do *not* do — hence the untested core of this report:**
1. **They name `spec Δ_W=[0,∞)` but never use the *type*.** They argue "the wedge flow is the free-field **boost**, therefore mass-blind" — a statement *about the free field*. The stronger, **theory-independent** statement is: `σ(Δ_W)=[0,∞)` *because `𝒜(W)` is Type III₁*, and the type is fixed by the net+vacuum, so the spectrum is mass-blind **even for interacting YM**, independently of BW's geometric identification. This is the actual no-go backbone and it is absent. (New §2.)
2. **They never run the anti-localization test on the modular route as a whole.** `AQFT-CAUSAL §4` notes QNEC is "wrong channel" for the delocalized `k→0` edge; none derives the general mechanism `K_W=2π∫x¹T₀₀ ⇒ modular is anti-localized off the sign cell`, nor states modular theory as the ninth costume of `TRANSPORT-CORE §4`'s anti-localization principle. (New §4.)
3. **They assert the BW deviation is circular but do not separate interaction from mass.** New §4 makes the precise "two theories, same `m`, different deviation" / "bounding `δK` reads the low-`k` stress tensor" argument, and kills the untested sub-levers (finite-region ball modular Hamiltonian, relative-entropy first law). (New §4.)

So: modular theory is a re-tread of the banked AQFT cluster **on the mass-blind-boost axis**; what is new is the **Type III₁ spectral theorem as a clean, theory-independent no-go**, plus the proof that modular **reproduces** rather than evades anti-localization. Both are negative — that is the honest result.

---

## §1. Modular structure, stated exactly (真定理 / standard)

**The net and the vacuum.** A Haag–Kastler net `O ↦ 𝒜(O)` of von Neumann algebras on a Hilbert space `ℋ`, with a unique Poincaré-invariant vacuum vector `Ω`. **Reeh–Schlieder**: `𝒜(O)Ω` is dense for every open `O` — equivalently `Ω` is **cyclic and separating** for each `𝒜(O)` (and for each wedge `W`).

**Tomita–Takesaki (真定理; Tomita 1967 / Takesaki 1970).** For `(ℳ,Ω)` with `Ω` cyclic+separating there is a unique **modular operator** `Δ_ℳ` (positive, self-adjoint, nonsingular) and **modular conjugation** `J_ℳ` (antiunitary),
```
Δ_ℳ = J_ℳ ∇_ℳ = ∇_ℳ J_ℳ,    J_ℳ Δ_ℳ J_ℳ = Δ_ℳ⁻¹,    Δ_ℳ Ω = Ω,  J_ℳ Ω = Ω,
```
such that the **modular automorphism group** `σ_t = Ad Δ_ℳ^{it}` preserves `ℳ` (`σ_t(ℳ)=ℳ`), and `J_ℳ ℳ J_ℳ = ℳ′` (the commutant). The **modular Hamiltonian** is `K_ℳ := −log Δ_ℳ`, so `σ_t = Ad e^{−itK}` is "time evolution under `K`" and `K Ω = 0`.

**KMS = thermal at `β=2π` (真定理).** The vacuum is a KMS state for the modular flow at inverse temperature `β=2π` (Araki; Haag–Hugenholtz–Winnink): for `A,B∈ℳ`, `ω(AB)=ω(B σ_{i2π}(A))`. This is the **algebraic** content — no measure, no positivity cone, just cyclicity + the modular condition. *(This is the honest sense in which the species is new: it is KMS/algebra, disjoint from the reflection-positivity/FKG/Γ₂ family.)*

**Bisognano–Wichmann (真定理 for free; conjectural-in-UV for interacting).** For the **Rindler wedge** `W={x:x¹>0}`, the modular flow is the **geometric boost** `Δ_W^{it} = U(Λ_W(−2πt))` and `J_W = CPT·R` (a `π`-rotation `R` about the `x¹`-axis composed with CPT). Equivalently `K_W = 2π K_boost = 2π ∫_W x¹ T₀₀(x) d³x`. For a general interacting theory, `K_W` is expected to equal the boost in the **UV** (leading area-law) with **subleading corrections** `δK_W = K_mod − K_boost ≠ 0` — this `δK_W` is the **BW deviation**. (Sewell; the free-field exact case is Bisognano–Wichmann 1975/76.)

**The type (the load-bearing fact).** Local algebras of a 4D QFT are **Type III₁** — established for the free field and for asymptotically-free models (Fredenhagen; Yngvason), and the **expected universal type** for 4D YM (conditional on the not-yet-constructed YM net, the G1 wall). Connes' classification fixes the **modular spectrum / `S`- and `T`-invariants**:
```
Type III₁:   S(𝒜) = spectrum of Δ = [0,∞),     T(𝒜) = ℝ   (modular flow has no periodicity),
```
i.e. `σ(Δ_W)=[0,∞)`, equivalently `σ(K_W)=−log[0,∞)=ℝ` — the task's "spectrum is the full line `ℝ`." There is **no smallest positive spectral value** (`0∈σ(Δ)` as a limit point), the modular flow is **ergodic** (no almost-periodic part; the only `Δ`-invariant vectors in the standard form are trivial), and `𝒜(W)` has **no normal invariant state and no trace**.

**The target being attacked (H1).** `Clay gap ⟺ Var_μ(a_{k→0}) ≤ 1/λ`, `(L,a)`-uniform, equivalently the atom `f₁²=|⟨0|trF²|G⟩|²>0`. In modular language the forge is: *use `(𝒜(W),Δ_W,J_W,Ω,KMS)` to bound `Var_ω(a_{k→0})` from above.* §2–§4 test this.

---

## §2. The Type III₁ gapless-spectrum obstruction (the CORE no-go)

This is the precise new statement. It is **stronger** than the banked "the free-field boost is mass-blind," because it does **not** rely on BW's geometric identification and holds for the **interacting** theory as well.

### §2.1 Theorem (Type III₁ ⇒ no modular gap ⇒ no mass bound) — 条件定理, conditional only on `𝒜(W)` being Type III₁

Let `𝒜(W)` be the wedge algebra with vector `Ω` (cyclic+separating), modular operator `Δ_W`, modular Hamiltonian `K_W=−logΔ_W`. Suppose `𝒜(W)` is Type III₁. Then:

**(a) `σ(Δ_W)=[0,∞)`, `σ(K_W)=ℝ`** (Connes' `S`-invariant of Type III₁) — the modular spectrum is the **full half-line / full line**, with `0` a non-isolated limit point and **no spectral gap**: `inf{λ>0 : λ∈σ(Δ_W)} = 0`, `inf |σ(K_W)| = 0`.

**(b) The modular flow `σ_t=AdΔ_W^{it}` is ergodic** — it has no almost-periodic / no point-spectrum component; `Δ_W` has no eigenvalue `≠1` with a vector in the physical standard form, and there is **no invariant normal state** (`𝒜(W)` is not semifinite, has no trace, no minimal projections).

**(c) `σ(Δ_W)` is mass-independent.** The type is determined by the **net + vacuum** (a UV/kinematic datum: Reeh–Schlieder + weak additivity + the continuum of local degrees of freedom), *not* by the mass gap. The free massless, free massive, and (expected) gapped interacting YM wedge algebras are **all Type III₁**, so `σ(Δ_W)=[0,∞)` for **all** of them.

**Consequence (this is the no-go).** Any quantity `F(Δ_W)` that depends **only on the spectrum / spectral type** — e.g. "the modular gap `inf spec Δ_W`", "the bottom of the modular spectrum", "the number of discrete modular eigenvalues", "the spectral dimension of the modular flow", "the modular `S`/`T` invariant" — takes the **same mass-independent value** on the gapless and gapped theories. In particular `inf{λ>0:λ∈σ(Δ_W)}=0` **identically**, so *a naive inequality "modular spectral gap `ε>0` ⇒ mass gap `m>0`" has a false premise*: the modular gap is `0` even when `m>0`. The modular **spectrum** is a mass constant and cannot lower-bound `m`.

### §2.2 The robustness dichotomy — why the no-go survives the absence of a constructed YM net

One might object: "Type III₁ for YM is conditional (no net yet, G1 wall)." The no-go is **robust regardless**, by a structural dichotomy.

**(Dichotomy — 条件定理, on standard QFT axioms).** For a 4D QFT local algebra `𝒜(W)` (Reeh–Schlieder, weak additivity, a Lorentz-covariant vacuum):
- **Branch III₁** (the expected/established one): `σ(Δ_W)=[0,∞)`, gapless, ergodic ⇒ by §2.1 **no modular spectrum lever**.
- **Branch not-III₁** (Type I or II): Type I ⇒ `𝒜(W)≅B(ℋ_W)` has **minimal projections** ⇒ the region carries finitely many degrees of freedom, contradicting the continuum (arbitrarily many independent local modes; the split property's infinite nuclearity). Type II ⇒ `𝒜(W)` has a faithful normal **semifinite trace** ⇒ the modular flow has a **periodic/discrete** spectrum (`T`-invariant `{0}`), and an invariant weight — **incompatible** with the `β=2π` KMS vacuum of a Lorentz-covariant field (whose modular flow must be geometric and, by BW, a boost with spectrum `ℝ`). So Type I/II is **excluded** for a genuine 4D continuum field.

**Therefore:** either `𝒜(W)` is Type III₁ (⇒ §2.1, spectrum blind) **or** the modular framework is **inapplicable** (the algebra is not that of a 4D continuum QFT). In neither branch is there a mass **lower** bound from modular spectrum. **The no-go does not need G1.**

### §2.3 What the spectrum no-go does *not* claim (honesty)

The theorem kills the **spectral-set** route: `σ(Δ_W)` and every functional of it alone. It does **not** by itself say "every functional of `Δ_W` is mass-blind" — the modular **spectral measure** of a vector (`⟨ξ,E(·)ξ⟩`) can in principle carry mass information. That residual is exactly the **entanglement-spectrum** route, and it is dealt with honestly in §4.3: it is (i) a genuine discriminator **only in the Gaussian/free regime** (banked `MODULAR-TAKESAKI`: `S` saturates ⟺ `m>0` for free fields), and (ii) **requires the correlators** (= H1) for interacting YM (Li–Haldane). So the precise output is a **two-pronged** statement: *spectral-set functionals are mass-blind (Theorem §2.1); spectral-measure functionals are free-only or circular (§4.3)*. Together: **no modular lever.**

---

## §3. The two acid tests (both run decisively)

### Acid test (1) — does it avoid measure-positivity *and* anti-localization?

This is the crux, and it has a **split answer** — the genuine value of the species and its kill.

- **Avoids measure-positivity — YES (genuine).** The dead positivity/curvature species (`Γ₂`, reflection positivity, FKG, Cheeger, Perron) all live on **positivity of a measure / curvature**. Modular theory lives on the **algebra**: `(ℳ,Ω)↦(Δ,J)`, KMS at `β=2π`. It uses **cyclicity + the modular condition**, not a positivity cone. So this is a genuinely **different species** on the measure-positivity axis — the first in the campaign that is *algebra, not geometry*. This is why it earns serious examination rather than instant dismissal.

- **Avoids anti-localization — NO (the kill).** The modular Hamiltonian is `K_W=2π∫_W x¹ T₀₀(x)d³x`. It is an integral of the **local stress-energy density**. The sign-bearing cell (`TRANSPORT-CORE §4`: typical field × `k→0` × `O(π)` phase) is carried by the **delocalized `k→0` transverse soft mode**: amplitude `∼1/√V` spread over the whole volume, energy `ω(k)∼k→0`. Its **local energy density `T₀₀→0`** on that mode (a soft collective mode of total energy `ω/V→0` per unit volume). Therefore `K_W` — and every modular functional built from `∫x·T₀₀`, including `δK_W` — is **suppressed** on the sign cell and **maximal on the UV / area-law cell** (where `T₀₀` is dominated by local high-momentum modes). **Modular structure is anti-localized to exactly the cell where the gap is *not* decided.** (§4 makes the general table.)

**Result: (1) PASS on measure-positivity, FAIL on anti-localization.** Modular theory is a new species that **reproduces** the campaign's universal support-mismatch wall (§4). It did not escape; it re-entered through `T₀₀`.

### Acid test (2) — free massless field (gapless, `f₁²=0`): no spurious gap?

- **Wedge / spectral route.** Free massless: `Δ_W^{it}=boost` exactly (BW), Type III₁, `σ(Δ_W)=[0,∞)`, `K_W=2πK_boost` with `K_boost` the boost; `Δ_W Ω=Ω` (`K_WΩ=0`). The spectrum is **gapless and continuous** — the route returns **no mass (`0`)**. `f₁²=0` (free Maxwell has positive-definite PF vacuum yet no atom; `PROOF-GAUSS-INVARIANT-PF`). **No spurious gap: PASS.** But — by §2.1(c) — the **same** spectrum holds for free massive and for (expected) gapped YM: the spectral route returns `0` for **all** masses. So the free-massless pass is **correct but non-discriminating** (vacuous on the mass axis).
- **Entanglement-entropy route (free/Gaussian).** Banked `MODULAR-TAKESAKI` (numpy, 1D free-fermion chain, `L≤800`): `S=(1/6)lnL+0.287` at `m=0` (log-divergent, **no saturation** ⇒ no gap), `S=0.7292` / `0.4887` (`L`-independent) at `m=0.1` / `0.4` (**saturates ⇒ gap**). So the **finite-region** spectral-**measure** functional *does* correctly fail the free massless field **non-vacuously** and tracks the free gap. **PASS — but Gaussian-only** (the interacting extension is Li–Haldane = circular, §4.3).

**Result: (2) PASS.** No spurious gap on the free massless field; the Gaussian entanglement diagnostic even discriminates massless/massive correctly. The decisive nuance (`TRANSPORT-CORE §3.2`'s "half-vacuous" lesson): the spectral route is `0` for **all** masses (mass-blind, §2.1), and the entropy route is a genuine discriminator **only on the free axis**. The free-field pass therefore does **not** touch the YM target — confirming, not refuting, the no-go.

---

## §4. Does modular structure evade anti-localization? (NO — it reproduces it)

### §4.1 The `K_W=2π∫x¹T₀₀` anti-localization mechanism

The modular Hamiltonian of a wedge is, by BW and its interacting corrections, a **weighted local stress-tensor integral**:
```
K_W = 2π ∫_W x¹ T₀₀(x) d³x  (+ UV-finite edge corrections),      δK_W = K_mod − K_boost.
```
Every **spectral-measure** modular functional (entanglement entropy `S_W`, relative entropy `S(ω‖ω_vac)`, the "modular temperature" response, QNEC `⟨T_kk⟩`) is a functional of `K_W`, hence of `∫x·T₀₀` — i.e. a **local energy-density integral**. The sign-bearing cell σ* (typical field × `k→0` × `O(π)` phase) is the **delocalized soft mode**: `amp∼1/√V`, `ω∼k→0`, `T₀₀|_{σ*}→0`. Thus:

| modular functional | large where | value on the sign cell σ* (typical × `k→0` × `O(π)`) |
|---|---|---|
| `K_W`, `S_W`, `S_rel` | UV / area-law (local `T₀₀` dominant) | `→0` (delocalized soft mode: `T₀₀→0`) |
| QNEC `⟨T_kk⟩` (AQFT-CAUSAL §4) | local / lightlike-near-UV | `→0` (wrong channel: not the `k→0` edge) |
| BW deviation `δK_W` | interacting / UV | controlled by low-`k` `⟨T T⟩` = mass data (circular, §4.2) |
| `σ(Δ_W)`, type, `S`/`T`-invariant | everywhere (mass-constant) | `[0,∞)` — **no mass content at all** (§2.1) |

**Every modular handle is largest where the gap is not decided and vanishes (or is constant) on the one cell where the sign lives.** This is `TRANSPORT-CORE §4`'s anti-localization principle, now for the **ninth** time: `‖Hess f‖²` (high `k`) → `Φ″` surplus (log-cut `F≈0`) → WSTC (high field) → **`K_W=∫x·T₀₀` (UV/area-law)**, all vanishing on typical-field × `k→0` × `O(π)`. Modular theory is **algebraic, not geometric** — but its generator is still the **local** boost `∫x·T₀₀`, so it inherits the same support mismatch. **The species does not evade anti-localization; it is anti-localization viewed through `T₀₀`.**

### §4.2 The BW-deviation circularity, made precise (interaction ≠ mass)

The only modular object that is **not** the bare geometric boost is `δK_W=K_mod−K_boost` (nonzero when the theory is interacting). The hope: "`δK_W≠0` ⟺ interacting YM ⟺ gapped ⇒ bound `m` from `‖δK_W‖`." Why it fails:
- **`δK_W` measures interaction, not mass.** `δK_W` is fixed by the **connected stress-tensor correlators in the wedge** `⟨T_{μν}(x)T_{ρσ}(y)⟩_c`. Two theories can have the **same mass `m`** and **different `δK_W`** (differing UV/interaction content at fixed pole), and conversely. So `δK_W` is not a clean mass functional — bounding it gives an **interaction** bound, which the pure-YM continuum has at all couplings, gapped or not.
- **Bounding `‖δK_W‖` reads the low-`k` stress tensor = the KL spectral measure = `m`.** Any quantitative bound on `δK_W` flows from the low-momentum `⟨TT⟩_c`, whose Källén–Lehmann representation is `∫ ρ(s)/(s−k²) ds` with `inf supp ρ = m`. So controlling `δK_W` at `k→0` **is** controlling the mass residue — i.e. `δK_W`-bound ⇒ KL-spectral-measure bound ⇒ H1. **Circular in content** (exactly the `MODULAR-TAKESAKI`/`AQFT-CAUSAL` wall, here localized to `δK_W`).

So even the *interacting-corrected* modular Hamiltonian is no escape: its deviation is an interaction functional, and the only route to a lower bound on it passes through the mass data it is supposed to produce.

### §4.3 The untested sub-levers, and why each dies (no lazy no-go)

To make the no-go complete rather than dismissive, the two genuinely-untested modular sub-levers are run and killed:

- **(L1) Finite-region (ball) modular Hamiltonian.** For a ball of radius `R` in a **CFT**, `K_ball=2π∫ (R²−r²)/(2R) T₀₀` (Casini–Huerta–Myers). A **massive/interacting** theory corrects it — and the correction **is** mass-sensitive (the finite-region entanglement spectrum sees the gap, for Gaussian systems). **Kill:** (i) the correction is still a **weighted local `T₀₀` integral** ⇒ the §4.1 anti-localization argument applies verbatim (the soft delocalized `k→0` mode has `T₀₀→0`); (ii) the ball-`S` mass-sensitivity **is** the entanglement entropy's, which for **interacting** 4D YM is the Li–Haldane conjecture (proven only Gaussian) = H1. **Dies exactly as the wedge route.**
- **(L2) Relative-entropy first law `δS(ω‖ω_vac)=δ⟨K_W⟩`.** The first law is an **identity**, not a bound: it relates the entropy change to the energy change, but to *use* it one must **know** `ω`'s energy (`=m`). And `S_rel` is **area-law/UV-dominated** (banked `MODULAR-TAKESAKI §4a`). **Kill:** circular (needs the state's mass) + wrong energy end (UV). It cannot produce a lower bound; it presupposes the answer.
- **(L3) Modular spectral measure / entanglement entropy, interacting.** The lone survivor of §2.3: a spectral-**measure** functional could in principle be mass-sensitive. **Kill:** established only for **Gaussian** systems (free fields); the interacting 4D YM version is Li–Haldane, which **is** the gap. (This is the precise residue; see exact stuck step.)

**Conclusion of §4.** No modular sub-lever evades anti-localization, and the only mass-sensitive one (entanglement) is Gaussian-only / circular. The modular species **reproduces** the anti-localization wall, and the honest output is the clean no-go of §2.

---

## §5. Scorecard (0–5), sharpest line, exact stuck step

### Five-condition grade: **1/5**

| # | condition | grade | one-line basis |
|---|---|---|---|
| 1 | **nonperturbative** | **PASS** | Tomita–Takesaki / KMS / BW are exact, all `β`; no loop, no `1/N`. Genuine algebra (disjoint from measure-positivity). |
| 2 | **sign-controlling LOWER bound `m_eff²(k→0)≥λ`** | **FAIL** | `σ(Δ_W)=[0,∞)` (Type III₁) is **mass-independent** (§2.1); no modular spectral gap; `K_W=boost` mass-blind; nothing gives `λ>0`. |
| 3 | **non-circular** | **FAIL** | Spectral route vacuously `0`; the only mass-sensitive functional (entanglement) needs the correlators (Li–Haldane = H1); BW deviation `δK_W` reads the low-`k` `⟨TT⟩` = KL measure = `m` (§4.2). |
| 4 | **non-SUSY + fails free massless field** | **PASS (vacuous on YM)** | No spurious gap on free massless (`S` log-divergent, `f₁²=0`, correctly `0`); Gaussian entropy even discriminates massless/massive. But spectral route is `0` for **all** masses ⇒ not a YM discriminator. Acid test 2 clean. |
| 5 | **`(L,a)`-uniform `m_eff²(k→0)≥λ`** | **FAIL** | `[0,∞)` is trivially uniform and gives `0` for every mass; no `λ`. Uniformity is vacuous. |

**Score = 1/5** (only #1 genuine; #4 passes the no-spurious-gap acid test but is vacuous as a YM discriminator; #2/#3/#5 fail). **Identical to the three banked modular files (`MODULAR-TAKESAKI` 1/5, `H1-T2-AQFT` 1/5, `AQFT-CAUSAL` 1/5)** — confirming, at the level of the *type* rather than the boost, that modular theory is a re-tread of the AQFT cluster. It is the **most principled** 1: it passes acid test 2 with the *correct* physics (a real free-field gap detector), and its kill is a **structural theorem of operator algebras** (§2.1), not a missing estimate.

### Sharpest line
*The modular species is the first that is **algebra, not geometry** — it genuinely escapes the measure-positivity family — and that is precisely why its failure is the **cleanest**: because the local algebras are **Type III₁**, the modular spectrum is `σ(Δ_W)=[0,∞)` (equivalently `σ(K_W)=ℝ`) **for every mass**, so the "modular gap" is `0` even when `m>0` — a spectral-set functional is a **mass constant**, and the only non-constant handles (`K_W=2π∫x¹T₀₀`, the BW deviation `δK_W`) are a **local stress-tensor integral** that vanishes on the **delocalized soft `k→0` sign cell**, and an **interaction functional** whose bound reads the low-`k` `⟨TT⟩` = the KL measure = `m`. The ninth costume of anti-localization, sewn in algebra: modular theory sees the **geometry of the vacuum**, never the **stiffness of its softest mode**, and when it tries, it must first know the mass it is trying to bound.*

### EXACT stuck step (缺失步)
```
target:   m_eff²(k→0) ≥ λ  uniform (L,a)   ⟺   Var_μ(a_{k→0}) ≤ 1/λ   =  H1.
needed:   a functional  Φ[ 𝒜(W), Δ_W, J_W, Ω ]  with
             Φ = 0 on the free MASSLESS vacuum (gapless),   Φ ≥ λ>0 on gapped YM,
             and Φ NOT reducible to  (i) the spectrum/type of Δ_W,
                                     (ii) a local  ∫ x·T₀₀  integral,
                                     (iii) the low-k stress-tensor correlator.
NONE exists on the existing modular algebra:
  (i)  spectrum:  σ(Δ_W)=[0,∞)  is Type III₁,   MASS-INDEPENDENT            → 0 for all m  (定理 §2.1)
  (ii) local:     K_W=2π∫x¹T₀₀  is anti-localized to the UV/area-law cell; on σ* (delocalized k→0,
                  T₀₀→0) it vanishes                                            → 0 on the sign cell  (§4.1)
  (iii) correlator: δK_W / entanglement = low-k ⟨TT⟩ = KL spectral measure, inf supp = m  → circular = H1  (§4.2)
missing:  a mass-sensitive modular functional that ALSO lives on the sign-bearing cell
          and does NOT read the correlators — i.e. an operator-algebraic coercivity on the
          delocalized k→0 mode.  This is H1 in KMS dress.  It does not exist.
```

### What would change the verdict (honest)
Not a sharper modular estimate — a change of **type**: (i) a proof that `𝒜(W)` is **not** Type III₁ (would mean the region has a trace / discrete modular spectrum — in tension with Reeh–Schlieder + the `β=2π` KMS vacuum; the dichotomy §2.2 says this is excluded for a 4D field, and if it held the modular framework itself would be inapplicable), or (ii) a **mass-sensitive** modular inequality on the `k→0` mode that avoids the local `T₀₀` integral (§4.1) and the correlators (§4.2) — i.e. an algebraic coercivity that does not read `⟨TT⟩`. Absent one of these, the 10th species is closed as a **clean no-go**: *modular theory is blind to the mass gap because Type III₁ makes the modular spectrum gapless independent of mass, and the only mass-sensitive modular data is circular.* That no-go is itself a precise, reusable result — the operator-algebraic twin of `TRANSPORT-CORE §4`'s anti-localization principle.

---

*Refs.* **Internal (banked, cited as ground, not re-derived):** **MODULAR-TAKESAKI** (free entanglement diagnostic `S~(1/6)lnL` massless vs saturation massive; wedge boost mass-blind; interacting ES = Li–Haldane = gap; relative entropy UV, 1/5); **H1-T2-AQFT** (AQFT = statement layer; modular kinematic mass-blind; mass is input; QNEC local/wrong-channel, 1/5); **AQFT-CAUSAL** (modular-mass-blindness lemma `Δ_W^{it}=boost`, `spec[0,∞)`, `β=2π`, same `m=0`/`m>0`; free massless net = counterexample; RS anti-clustering; split = UV; clustering rate = H1, 1/5); **TRANSPORT-CORE §4** (anti-localization principle: positivity largest where gap not, vanishes on typical-field × `k→0` × `O(π)` cell — the sign cell σ*); PROOF-SKELETON §0 (量纲 `G″_a=a²m²`), §5 (single species = variance upper bound, UV=L3 / IR=H1); REDUCTION-AUDIT (H1 = susceptibility, `f₁²` atom is strictly stronger = `H*`); PROOF-GAUSS-INVARIANT-PF (free Maxwell PF vacuum yet `f₁²=0`); CLAUDE.md §1–2 (造新工具, 诚实铁律, free-field acid test). **External (standard):** Tomita 1967; Takesaki, *Tomita's theory of modular Hilbert algebras*, LNM **128** (1970); Connes, *Une classification des facteurs de type III*, Ann. Sci. ENS **6** (1973) 133 (S/T invariants, III₁ `S=[0,∞)`, `T=ℝ`); Connes, *Almost periodic states and factors of type III₁*, J. Funct. Anal. **16** (1974) 415; Bisognano–Wichmann, J. Math. Phys. **16** (1975) 985 / **17** (1976) 303; Araki, PTP **32** (1964) 844; Haag–Hugenholtz–Winnink, CMP **5** (1967) 215; Reeh–Schlieder, Nuovo Cim. **22** (1961) 1051; Fredenhagen, *On the modular structure of local algebras of observables*, CMP **97** (1985) 79; Yngvason, *The type of local algebras*, Lett. Math. Phys. **4** (1980) 363; Sewell, *Quantum fields on manifolds: PCT and gravitationally induced thermal states*, Ann. Phys. **141** (1982) 201; Casini–Huerta–Myers, JHEP **05** (2011) 036 (ball modular Hamiltonian); Li–Haldane, PRL **101** (2008) 010504 (entanglement-spectrum gap ⟺ Hamiltonian gap, Gaussian; conjectural interacting); Unruh, Phys. Rev. D **14** (1976) 870.
