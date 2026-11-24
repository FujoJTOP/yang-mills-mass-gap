# TRANSPORT-SIZE-FIX — flip condition #2 of `E2-FLIP-TRANSPORT §5`: can the instanton SIZE MODULUS in ℝ⁴ be **fixed** by a mechanism derived from the action? Verdict: **NOT ACHIEVED** — a stationary point exists but only from one-loop **measure × running**, it sits at `O(1)/Λ` (not `≪1/Λ`) at `α_s ≈ 0.5` (one-loop uncontrolled), and the IĪ "molecule/liquid" stabilization is a **FIT** (self-consistent `{n, ρ̄}`), **circular in `Λ`**; the acid test both **passes** and, in passing, **proves the circularity**.

> 2026-11-24 · Direct attempt on flip #2 of `E2-FLIP-TRANSPORT §5` (and the "machine that does not exist", `CLAUDE.md §5`): *construct or refute a size-fixing mechanism for the instanton modulus in ℝ⁴ that pins `ρ* ≪ 1/Λ` from the action, not fitted.* We write the effective action for the collective coordinate `ρ`, add the IĪ interaction both as a **mean-field liquid** (setting A) and as an **isolated molecule** (setting B), locate the minimum numerically for SU(3)/SU(2), run the **free-field acid test**, and separate **theorem from fit**.
>
> **Labels (`CLAUDE.md §2`):** **[T]** 真定理 · **[C]** 条件定理 · **[NG]** no-go *for a named class* · **[O]** 精确缺失步 · **[J]** 判断 · **[V-here]** fresh computation this session (`research/unconditional/w7/transport_sizefix.py`). External facts are **[R]** (as-reported; not re-verified here). **No QED on anything unproven. This note does NOT prove Clay and does NOT prove no mechanism can exist (`∀class ≠ ∀possible`).**

---

## Verdict up front

> ### **NOT ACHIEVED (flip #2 stays closed).** The size effective action **does** have a stationary point — but **only** from the one-loop competition between the **measure** `dρ/ρ⁵` (favours small `ρ`) and the **running action** `8π²/g²(ρ)=b\ln(1/ρΛ)` (favours large `ρ`). It sits at
> ```
>    ρ*  =  e^{−2N/(b−5)} / Λ  =  O(1)/Λ      (SU(3): 0.368/Λ ;  SU(2): 0.180/Λ),
> ```
> **at strong coupling** `α_s(ρ*) = 2π/(b\ln(1/ρ*)) ≈ 0.57 (SU(3)) / 0.50 (SU(2))` — i.e. exactly where the one-loop running used to *find* it has broken down. It is **not** at a parametrically small `ρ* ≪ 1/Λ`, and its **location is `c/Λ`**, so it *imports the very scale it is supposed to generate*.
>
> **The IĪ interaction does not help and does not produce a competing small-`ρ` minimum.** It is **attractive** (`ΔΓ_int<0`), **bounded**, and its own extremum is at `ρ ∼ 1/e·(1/Λ)`; it pulls the size **toward `≈0.37/Λ`** (LARGER for SU(2), coincident for SU(3)) — **never smaller**. At the physical packing it is a **fit-level self-consistency**, and at `O(1)` packing the integral runs to `R ∼ ρ`, where the large-`R` potential `V ∼ 1/R⁴` is **invalid**.
>
> **The acid test PASSES and proves the circularity in the same stroke:** with **no running** (fixed `g`, no `Λ`) the exponent is monotone (`Γ_free = −5\ln ρ`), there is **no minimum**, the modulus is conformal. Hence **the minimum's very existence IS the running IS `Λ`** — you cannot use "the size is fixed at `0.37/Λ`" to certify a gap of size `Λ` without already having `Λ`.
>
> **One line:** the modulus is **not** removed in ℝ⁴; the apparent "0.33 fm stabilization" is the one-loop peak `ρ* ≃ 0.37/Λ` evaluated with `Λ_QCD` — a **picture/fit**, not a theorem — and converting it into a proof requires control of the `ρ`-effective action at `α_s ≈ 0.5` (the multi-instanton liquid at packing `O(1)`), which is the **same crossover wall (`H1`)**.

---

## §0 The target, stated so it can be refuted

Flip #2 (`E2-FLIP-TRANSPORT §5`): *exhibit, for pure SU(N) YM on ℝ⁴, a **provable** size-fixing mechanism that removes the instanton modulus (not on `S¹`, not via a deformation) and pins `S_inst` bounded across the crossover.* Concretely, an effective action `V_eff(ρ)` — derived from the QCD action, not fitted — whose minimum is at a **finite, parametrically small** `ρ* ≪ 1/Λ`, with a **stable** (2nd-order) minimum whose **location does not presuppose `Λ`**.

Two things must be separated throughout:
- **"A minimum exists"** (a calculus fact about whatever `V_eff` you write) — call it existence.
- **"The minimum is the physical, controlled size, fixed independently of `Λ`"** — call it control + non-circularity.
The task asks for *both*; we will find existence **but not** control + non-circularity.

---

## §1 The effective action for the size collective coordinate `ρ`

### 1.1 Classical + one-loop (measure, running, determinant) `[T-formula / V-here]`

For one instanton in pure SU(N) YM, the one-loop measure is the standard 't Hooft–Bernard result `[R]`
```
   dμ/dρ d⁴z  ∝  (1/ρ⁵) · [8π²/g²(ρ)]^{2N} · e^{−8π²/g²(ρ)} · d⁴z ,        8π²/g²(ρ) = b·ln(1/(ρΛ)) ,   b = 11N/3 ,
```
where `1/ρ⁵` is the size measure (dilatation), `[8π²/g²]^{2N}` is the `4N`-zero-mode/determinant prefactor, and `e^{−8π²/g²}` is the saddle weight. Writing `L ≡ ln(1/ρ)` and dropping `ρ`-independent constants, the integrand is `exp(−Γ_1ℓ(L))` with the **one-loop effective action**
```
   Γ_1ℓ(L)  =  (b − 5)·L  −  2N·ln(b L) ,                       [T-formula]
   Γ_1ℓ'(L) =  (b − 5) − 2N/L ,     Γ_1ℓ''(L) = 2N/L²  > 0 .
```
This is the **only intrinsic** object in the problem. Its unique stationary point is `L* = 2N/(b−5)`, i.e.
```
   ρ*  =  e^{−2N/(b−5)} / Λ ,        8π²/g²(ρ*) = b L* = 2Nb/(b−5) ,     α_s(ρ*) = 2π/(8π²/g²) .
```
Three independent structural readings, all exact:
- **The measure is the ONLY term favouring small `ρ`** (`+5L = 5\ln(1/ρ)` grows as `ρ→0`).
- **Running and (below) the IĪ attraction both favour large `ρ`.** The minimum is measure-vs-running.
- **`ρ* ∝ 1/Λ`** — its *location is set by `Λ`*.

### 1.2 The IĪ interaction `[R]` (ratio ansatz / streamline)

The interaction of an instanton (`ρ`) and an anti-instanton (`ρ̄`) at separation `R` is (Callan–Dashen–Gross 1978; Schäfer–Shuryak RMP 1998 §III; Verbaarschot 1991 streamline; Diakonov–Petrov) **attractive** with the large-`R` **dipole–dipole** asymptotic `[R]`
```
   ΔS_{IĪ}(R;ρ,ρ̄)  =  − (32π²/g²) · w · (ρ² ρ̄² / R⁴)  +  O(ρ⁶/R⁶) ,     R ≫ ρ, ρ̄ ,
```
with `w` an order-unity orientation/colour factor. Two **limit facts** matter more than the coefficient and are coefficient-independent:
- **`R→0`:** `I` and `Ī` are topologically trivial together and can annihilate to the vacuum, so `ΔS → −2S₀ = −16π²/g²` (the **molecule** limit).
- **`R→∞`:** `ΔS → 0` (decoupling). Hence `ΔS` is **bounded below**, **attractive** (`ΔS<0`), and **saturates**.
*(The task brief sketched `V ∝ (32π²/g²)(ρρ̄/R²)`; the standard large-`R` form is the dipole–dipole `(ρ²ρ̄²/R⁴)` above — the `4D` field strength falls as `ρ²/r⁴`. Every conclusion below is coefficient-independent; only the sign and boundedness are used.)*

### 1.3 Setting A — mean-field liquid `[V-here]`

A probe instanton of size `ρ` in a frozen background of density `n` (size `ρ̄`) feels
```
   ΔΓ_int(ρ) = n ∫_{R_c}^{R_max} d⁴R  ΔS_{IĪ}(R;ρ,ρ̄) ,   d⁴R = 2π² R³ dR ,   R_c = max(ρ,ρ̄) ,   R_max = n^{−1/4} ,
```
where `R_max` is the inter-particle spacing (screening). Using `∫ dR/R = ln`, and `32π²/g² = 4·8π²/g²`
```
   ΔΓ_int(ρ)  =  − n · (32π²/g²(ρ)) · w · ρ² ρ̄² · 2π² · ln(R_max/R_c)   <  0     (attractive, bounded).
```
Two input regimes (Λ=1 units; `1 fm⁻¹ ≈ Λ_QCD ≈ 200 MeV`, so `1/Λ ≈ 1 fm`):

| regime | `n` (Λ⁴) | `ρ̄` (1/Λ) | packing `n ρ̄⁴` | source |
|---|---|---|---|---|
| **physical ILM** | `≈ 1` | `≈ 0.33` | `≈ 0.012` (≈1%) | Schäfer–Shuryak, Diakonov–Petrov `[R]` |
| **dense (free one-loop peak)** | `≈ 50` | `≈ 0.37` | `≈ O(1)` | `E2-FLIP-TRANSPORT §4.1` |

*(Note the honest nuance: the **physical** fitted ILM is **dilute** (`≈1%`), while the **free one-loop** peak implies a density `n ∼ ρ*^{−4} ∼ O(10–50)/Λ⁴`, packing `O(1)`. The interactions are what dilute the free runaway to the physical point — but only within the self-consistent fit; see §3. Either way the conclusion is the same.)*

### 1.4 Setting B — isolated IĪ molecule `[V-here]`

One member's size action with the partner **fixed** at separation `R` (no `d⁴R` integral — `V` enters directly):
```
   Γ_mol(ρ; R)  =  Γ_1ℓ(ρ)  +  ΔS_{IĪ}(R;ρ,ρ̄) ,     ΔS = −(4·8π²/g²(ρ))·ρ²ρ̄²/R⁴ .
```

---

## §2 The numerical minimum (or its absence) `[V-here]`

Script `research/unconditional/w7/transport_sizefix.py`, `py -3.12`, single process. `Λ=1`.

### 2.1 The intrinsic (one-loop) minimum — Table 1

| `N` | `b=11N/3` | `L*=2N/(b−5)` | **`ρ*`** | `8π²/g²(ρ*)` | **`α_s(ρ*)`** | `Γ''(L*)` | `σ_L` → size factor | 't Hooft `g²N` |
|---|---|---|---|---|---|---|---|---|
| SU(2) | 7.333 | 1.714 | **0.180/Λ** | 12.57 | **0.500** | 1.361 | 0.857 → ×2.4 | 13.6 |
| SU(3) | 11.00 | 1.000 | **0.368/Λ** | 11.00 | **0.571** | 6.000 | 0.408 → ×1.5 | 21.5 |

Readings: the minimum is **real** and **2nd-order** (`Γ''>0`), localising `ρ` to a factor `≈1.5–2.4`. But (i) `ρ* = O(1)/Λ`, **not** `≪1/Λ`; (ii) `α_s(ρ*) ≈ 0.5` — **strong**; (iii) the 't Hooft coupling at the peak is **large** (`g²N ≈ 14–22`), so even at large `N` the peak is at strong coupling (the instanton gas is dilute in `N`, but the *size* is fixed at strong 't Hooft coupling).

### 2.2 Adding the IĪ interaction — Tables 2 (mean-field) and 3 (molecule)

**Mean-field (Setting A):**
| `N` | regime | `ΔΓ_int` extremum at | its value | one-loop `ρ*` | **total `ρ*`** | shift |
|---|---|---|---|---|---|---|
| SU(2) | physical | `0.368` | −8.53 | 0.180 | 0.337 | **×1.87 (LARGER)** |
| SU(2) | dense | `0.370` | −8.76 | 0.180 | 0.370 | ×2.06 (LARGER) |
| SU(3) | physical | `0.368` | −12.80 | 0.368 | 0.368 | ×1.00 (coincides) |
| SU(3) | dense | `0.370` | −13.15 | 0.368 | 0.370 | ×1.01 (LARGER) |

**Isolated molecule (Setting B, `ρ̄=0.33`):**
| `R` (1/Λ) | `V(R,ρ̄)` | `dV/dρ` at `ρ̄` | interior min |
|---|---|---|---|
| 1.0 | −0.579 | −3.51 (attractive) | `ρ=0.405` (LARGER) |
| 2.0 | −0.036 | −0.22 (attractive) | `ρ=0.370` (LARGER) |
| 4.0 | −0.002 | −0.01 (attractive) | `ρ=0.368` (LARGER) |

**The decisive point:** `ΔΓ_int < 0` and **bounded**; its own extremum is at `ρ ∼ 1/e·(1/Λ) = O(1)/Λ`. Because it is **attractive and grows with `ρ`**, it pulls the size **UP toward `≈0.37/Λ`** (a big shift for SU(2) whose one-loop `ρ*=0.18`; coincident for SU(3) whose `ρ*=0.37`, since `2N/(b−5)=1` there). **In no regime does the interaction create or move a minimum to `ρ ≪ 1/Λ`.** The "molecule binding" lowers the pair action by **enlarging** the size until `R ∼ ρ` (the validity edge), not by pinning `ρ` small.

> **Why "attractive ⇒ larger":** the only force toward *small* `ρ` is the **measure** `1/ρ⁵`. Running **and** IĪ attraction both lower the action for larger `ρ` (a bigger instanton has a bigger dipole moment, stronger attraction, and a smaller `8π²/g²(ρ)`). Attraction can only ever **add** to the "large-`ρ`" side of the balance.

---

## §3 Theorem vs fit — the honesty analysis (the heart of the note)

We separate four claims that are routinely conflated in the ILM literature.

**(a) Existence of the stationary point — `[T]`, but vacuous as a *physical* fixing.**
The statement *"`Γ_1ℓ(L) = (b−5)L − 2N\ln(bL)` has a strict minimum at `ρ*=e^{−2N/(b−5)}/Λ`"* is **true calculus** `[T]` about the one-loop integrand. But the integrand is a **truncation** (one-loop). At its own minimum, `α_s ≈ 0.5`, so the truncated formula is **evaluated where it fails**. Calling the minimum "the physical instanton size" is therefore an **extrapolation to strong coupling**, not a theorem: `[O]`.

**(b) The location is `c/Λ` — `[T-tautology] ⇒ circular.**
`ρ* = e^{−2N/(b−5)}/Λ` is *one pure number times the only scale in the problem*. It does **not** fix the size relative to anything else; it says "the instanton size is `≈0.37` of the confinement scale." To use this to *produce* a mass gap of order `Λ` is to **assume `Λ`** — the scale the gap should generate. This is the same one-scale circularity banked as `E2-FLIP-NOGO §3 Theorem S` and `HIGGS-DIM-TRANS` ("transmutation fixes `Λ`, not the pole"). The minimum is **not at `ρ* ≪ 1/Λ`** — it is at `O(1)/Λ`.

**(c) The IĪ "stabilization at 0.33 fm" is a FIT, not a theorem `[R]/[J]`.**
- The inputs `n ≈ 1 fm⁻⁴`, `ρ̄ ≈ 0.33 fm` are **fitted** to lattice / phenomenology (Schäfer–Shuryak; Diakonov–Petrov variational ILM). The "0.33 fm" emerges from a **joint self-consistency** (stationarity of a model free energy in `{n, ρ̄}`), which is a **variational fit**, not a stability bound.
- At the relevant density the interaction is **not controlled**: the large-`R` form `V ∼ 1/R⁴` is valid only for `R ≫ ρ`, but at `O(1)` packing the dominant neighbours are at `R ∼ ρ` — the **validity edge**. So even the *correction* that would shift `ρ*` is computed where its own expansion fails.
- The interaction's sign is **destabilizing toward large `ρ`** (§2.2): it cannot be the agent that *fixes* a small size; at best it enlarges it.

**(d) The "0.33 fm = one-loop `ρ*`" is a numerical coincidence of the one-loop formula `[J]`.**
The physical ILM size `0.33 fm` and the one-loop peak `ρ* = 0.368/Λ_QCD ≈ 0.37 fm` agree to `≈10%`. This is *consistency*, not derivation: the one-loop formula happens to land on the strong-coupling scale where the real physics sits. Evaluating it with `Λ_QCD` and calling the agreement a "mechanism" is the postulate, not the conclusion.

**Net honesty verdict.** Existence: `[T]`. Controlled physical size-fixing: **`[O]`** (strong coupling at the min). Non-circular / `ρ* ≪ 1/Λ`: **fails**. IĪ molecule/liquid as the fixing agent: **fit `[R]/[J]`, sign-destabilizing, validity-limited.** **Flip #2 NOT achieved.**

---

## §4 The acid test (mandatory) `[V-here]`

Free / scale-invariant limit: **no running** (`g` fixed), **no `Λ`**. Then `8π²/g² = const` and the determinant prefactor is `const`, so the whole exponent is the pure measure
```
   Γ_free(L) = 5 L  =  −5\ln ρ ,      Γ_free'(L) = 5 > 0  everywhere ⇒ NO stationary point,
   dμ ∝ ρ^{−5} dρ   is scale-covariant:  ρ → λρ  multiplies it by  λ^{−4} .
```
Numerically (`transport_sizefix.py` Part C): no minimum, no preferred size, the size integral `∫ρ^{−5}dρ` is UV-divergent (dominated by the cutoff, `m_inst=0`). **ACID PASSED**: the mechanism **vanishes** on the free massless field — so it is **not** a gauge/regularization artifact.

> **The load-bearing consequence.** Passing the acid test and being **circular** are the **same fact** here: the *only* ingredient that created the PART A minimum was the **running**, i.e. **`Λ`**. Remove `Λ` → remove the minimum. Therefore *"there is a minimum at `c/Λ`"* **is** *"`Λ` is present."* A mechanism whose entire content is "the size is a fixed fraction of the confinement scale because the confinement scale is the only scale" cannot certify that the confinement scale is nonzero. The acid test is the **diagnosis**, not a passing grade for the mechanism.

---

## §5 The precise missing step, and the flip condition

The one statement whose proof would turn this from a fit into a theorem (equivalently flip #2):
```
   [O]  ∃ a controlled effective action  V_eff(ρ)  for the instanton size in pure SU(N) YM on ℝ⁴,
        valid (under analytic control, not one-loop extrapolation) on the range of ρ where its
        minimum sits, with minimum at a finite ρ* and a stability bound  V_eff''(ρ*) ≥ c > 0 ,
        uniform in (L,a),  WITHOUT importing Λ as an input.
```
By §1–§4 this is **not** a one-loop computation: the minimum sits at `α_s ≈ 0.5`, so `[O]` = **control of the multi-instanton measure / liquid at packing `O(1)` and `α_s = O(1)`** — the interacting instanton-liquid closure at the physical scale. That is **verbatim the crossover RG / `H1`** (`TRANSPORT-MACHINE Stage 3`, `E2-FLIP-TRANSPORT §4`). It is the same wall, reached through the size modulus.

**Why the IĪ route cannot shortcut it.** The IĪ interaction is `(i)` **attractive** (pushes to large `ρ`), `(ii)` only known in the `R ≫ ρ` asymptote (fails at `R ∼ ρ`, the packing edge), and `(iii)` an **ansatz** (ratio/streamline), not the path integral. So even a "molecule bound" would be a model statement, and its sign works against a small-`ρ` fixing. It cannot be the flip.

**Flip condition (stated so it can be refuted).** Flip #2 → **ACHIEVED** iff one exhibits a **provable** `V_eff(ρ)` for pure SU(N) YM on ℝ⁴ with a 2nd-order minimum at finite `ρ*`, **controlled** there (not a one-loop extrapolation to `α_s≈0.5`) and **non-circular** (the minimum and its scale derived from the dynamics, not written as `c/Λ`). None is known; the closest object, `ρ* = e^{−2N/(b−5)}/Λ`, fails control (strong coupling) and non-circularity (`c/Λ`).

---

## §6 Honest boundary — what this document does and does not say

- **Does NOT prove Clay; does NOT prove no mechanism can exist** (`∀class ≠ ∀possible`). It is a **`[V-here]` + `[T]/[O]/[J]` verdict on flip condition #2**: a genuine stationary point exists (`[T]`) but is uncontrolled (`[O]`) and circular (`[J]`); the IĪ "stabilization" is a fit (`[R]/[J]`), **not** a theorem, and is sign-destabilizing.
- **The positive facts are real.** (1) `Γ_1ℓ(ρ)` **does** have a strict minimum at `ρ* = e^{−2N/(b−5)}/Λ` `[T]`, localising the size to a factor `≈1.5–2.4`; this is *why* the one-loop size prediction is as good as it is. (2) The IĪ interaction is correctly attractive and bounded `[R]`; (3) the acid test passes cleanly `[V-here]`. Under-claiming these would be as dishonest as over-claiming the negative.
- **The `[V-here]` computation is a structural witness, labelled as such.** The one-loop minimum and the interaction scalings are evaluated with standard one-loop / ILM inputs (`[R]` CDG, Schäfer–Shuryak, Diakonov–Petrov); the numbers locate the peak, the coupling, and the interaction shift. **No SU(N) gap number is claimed.**
- **Relation to banked verdicts.** This sharpens `E2-FLIP-TRANSPORT §4.1` (the one-loop peak) into flip #2's exact question: the peak is `O(1)/Λ` at `α_s≈0.5`, the IĪ mechanism is a fit and pushes the *wrong* way, and the acid test **is** the circularity. It confirms `TRANSPORT-MACHINE Stage 3 = H1` in the size-modulus coordinate, and the one-scale circularity of `E2-FLIP-NOGO §3` / `HIGGS-DIM-TRANS`.
- **No fabricated mechanism, no QED.** Size-fixing in ℝ⁴ is **not** made rigorous here; it is shown to be **as hard as the gap** (control at `α_s=O(1)`) and **circular** in its only controlled-looking form.

---

## §7 References

**Internal (absolute paths):**
- `D:\Dev\math\CLAUDE.md` (§1 build-don't-stop + acid test; §2 labels; §5 wall = "the machine that does not exist"); `D:\Dev\math\docs\BRIEF-FOR-MODEL.md` (§4 "transport ≡ confinement").
- `D:\Dev\math\research\unconditional\E2-FLIP-TRANSPORT.md` (§4.1 one-loop size distribution + peak `ρ*`; §5 flip condition #2 — the exact item attacked here); `TRANSPORT-MACHINE.md` (4-stage machine; Stage 3 = `H1`); `E2-FLIP-NOGO.md` (§3 Theorem S, one-scale circularity); `W2-transmuted-gap.md` (gas needs diluteness + IR cutoff; ℝ⁴ fails); `HIGGS-DIM-TRANS.md` (transmutation fixes `Λ`, not the pole); `SMALLCIRCLE-ADIA.md` (size fixed only on `S¹`/deformation).
- `D:\Dev\math\research\unconditional\w7\transport_sizefix.py` (**the `[V-here]` witness**: Part A one-loop minimum + `α_s(ρ*)`; Part B mean-field interaction; Part B2 isolated molecule; Part C acid test). Run `py -3.12 transport_sizefix.py`.

**External — `[R]` standard / as-reported (not re-verified this session):**
- Belavin–Polyakov–Schwartz–Tyupkin, *Phys. Lett. B* **59** (1975) 85 (instantons, **size modulus**, classical `8π²/g²`).
- 't Hooft, *Phys. Rev. D* **14** (1976) 3432 (instanton measure, `[8π²/g²]^{2N}` prefactor); Bernard, *Phys. Rev. D* **19** (1979) 3013 (one-instanton measure, `1/ρ⁵`).
- Callan–Dashen–Gross, "The structure of the gauge theory vacuum," *Phys. Rev. D* **17** (1978) 2717 (IĪ interaction, large-`R` asymptote).
- Schäfer–Shuryak, "Instantons in QCD," *Rev. Mod. Phys.* **70** (1998) 323 (§III IĪ interaction / ratio ansatz; instanton-liquid size `ρ̄≈0.33 fm`, `n≈1 fm⁻⁴`; dilute-liquid packing).
- Shuryak, *Nucl. Phys. B* **203** (1982) 93 / 116 / 140 (instanton liquid); Diakonov–Petrov, *Nucl. Phys. B* **245** (1984) 259; Diakonov–Petrov in *At the Frontier of Particle Physics* (World Sci., 2000) (variational ILM: size from `{n,ρ̄}` self-consistency = a fit).
- Verbaarschot, *Nucl. Phys. B* **362** (1991) 33 (streamline / ratio-ansatz IĪ configuration).
- Coleman, *Aspects of Symmetry* (1985) Ch. 7 (tunnelling weight `e^{−S}`, rigid saddle); Kotecký–Preiss, *CMP* **103** (1986) 19 (cluster expansion / diluteness criterion — what the `O(1)` packing lacks).

---

*Status: **FLIP #2 NOT ACHIEVED.** A `ρ`-stationary point exists `[T]` but only from one-loop measure-vs-running, at `ρ* = e^{−2N/(b−5)}/Λ = O(1)/Λ` (SU(3): `0.368/Λ`, `α_s≈0.57`; SU(2): `0.180/Λ`, `α_s≈0.50`) — uncontrolled `[O]` and circular `[J]`. IĪ mean-field `[V-here]` and isolated molecule `[V-here]` are attractive, bounded, extremum at `∼1/e·(1/Λ)`, pulling the size LARGER — never to `ρ ≪ 1/Λ`; the "0.33 fm" stabilization is a self-consistent fit `[R]/[J]`. Acid test PASSED `[V-here]` and thereby PROVES the circularity (no running ⇒ no minimum). Exact missing step `[O]` = control of the `ρ`-effective action at `α_s ≈ 0.5` / packing `O(1)` = the crossover RG / `H1`. NOT a proof of Clay; NOT a global impossibility (`∀class ≠ ∀possible`). **Clay stays open.***
