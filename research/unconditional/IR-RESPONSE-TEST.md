# IR-RESPONSE-TEST — does the IR-non-analytic part of the connected 3rd-order nonlinear response `𝔑_IR^(n)` distinguish gapped from gapless 4D pure Yang–Mills?

> 2026-11-24 · **Decisive, honest evaluation of a candidate "6th species"** (charter §1). Candidate `𝔑_IR^(n)` = the IR-non-analytic part of the connected n-th-order causal nonlinear response (n=3), smeared `𝔑 = lim_{Λ→0} Λ^{−σ}∫d³u K(u) χ^{(3)}_{nonan}(Λu)`. Claimed property: **gapped ⇒ analytic at ω=0 ⇒ `𝔑=0`; gapless ⇒ non-analytic ⇒ `𝔑≠0`**, offered as a non-positivity, nonabelian, IR-support-sensitive gapped/gapless discriminator.
>
> Witness `research/unconditional/w7/ir_response_test.py` (ran clean, `py -3.12`, single job, ~9 s, numpy only), results `…/ir_response_test_results.json`. **No QED, no fake closure, no false hope.** This does **not** decide Clay; it decides whether `𝔑` is a usable lever.

**Labels (charter §2).** **[T]** true theorem · **[C]** conditional · **[O]** exact missing step · **[NG]** no-go for a named class · **[V-here]** computed this round · **[R]** external, as-reported/spot-checked.

---

## VERDICT (up front)

> ### The functional is **NOT a live handle**. It is **(a) a REFORMULATION of the gap** on the linear horn, **free-vacuous on the nonlinear horn** (Gaussianity kills it in the *mislabeling* direction), and **(b) non-discriminating** if derived from the trace anomaly. It is **not (c)**.
>
> **The dichotomy is genuine and collapses against the tool — in three independent ways:**
> 1. **Linear response** — "analytic at ω=0 ⟺ gap" is **correct** `[T]` but it *is* the Källén–Lehmann spectral representation. Deciding it requires `min supp ν`, i.e. the gap. **Reformulation ⇒ circular as a tool.**
> 2. **Nonlinear response (n=3)** — a free field is **Gaussian**, so **every connected n≥3 response is identically 0** for *all* masses `[T]`. Hence `𝔑^(3)(free massless) = 0 = ` the *gapped* value. The detector **mislabels the canonical gapless theory as gapped** and **fails the free-field acid test** (`[T]` + `[V-here]` Parts 1–2).
> 3. **Trace-anomaly pinning** — the anomaly / dilaton sum rule fixes only the **convergent** moment `∫ρ(s)/s² ds`, which is finite and positive for **both** a gapped measure and the gapless witness `dν=s²e^{−s/Λ}ds` ⇒ **blind** `[V-here]` Part 4. So the computable (anomaly-derived) `𝔑` is **non-discriminating**; the discriminating (full-IR) `𝔑` is **circular**. It inherits the blindness.
>
> **Precise resolution of the dichotomy:** no single definition of `𝔑` is simultaneously *(i) discriminating*, *(ii) non-circular/computable*, and *(iii) correct on free-gapless*. Linear = (i)+(iii) but not (ii). Nonlinear = (ii) but not (i)/(iii) on free-massless. Anomaly = (ii) but not (i). A genuine "6th species" must extract gap information from a channel that is **neither** the connected-response non-analyticity (Gaussian-blind) **nor** the convergent anomaly moment (gap-blind). **The proposal uses exactly those two channels, so it is not the 6th species.**

---

## §1 · The dichotomy: linear vs nonlinear

### 1.1 Linear horn — faithful but *= Källén–Lehmann* (reformulation)

`[T]` Write `G_R(ω)=∫_0^∞ dν(s)/(ω²−(s+i0)²)`. Then `G_R` is **analytic in a neighbourhood of `ω=0` iff `0∉supp ν`, i.e. iff the theory is gapped**. Equivalently (Part 3, verified numerically) define the KL moments
```
M_n := ∫_0^∞ ρ(s) s^{−2n−2} ds ,     G_R(ω) = −Σ_{n≥0} ω^{2n} M_n  (|ω|<m₁).
```
`G_R` is analytic at `ω=0` **⟺ every `M_n` finite** `[T]`. Gapped (`ρ=δ(s−m₁)`): `M_n=m₁^{−2n−2}<∞` ∀n. Gapless witness (`ρ(s)=s²e^{−s/Λ}`): `M_0=Λ` finite but `M_1=∫e^{−s/Λ}s^{−2}ds=+∞`.

**Numerical confirmation `[V-here]` (Part 3):** `M_1` gapped `= 4.165`; gapless witness `M_1 = 94.97, 1.01×10⁴, 1.5×10⁷` at IR cutoff `ε = 10⁻²,10⁻⁴,10⁻⁶` — divergent `~1/ε`.

**Verdict on the linear horn `[T]`.** The statement is true, but it is **the spectral measure restated**. To *evaluate* "is `G_R` analytic at 0" you must know whether all `M_n` are finite, i.e. whether spectral weight reaches 0, i.e. **the gap**. The functional is a **reformulation**, and using it as a tool is **circular**.

### 1.2 Nonlinear horn — Gaussian ⇒ identically zero ⇒ acid-test failure

`[T]` For a free scalar of **any** mass, `W=log Z` is quadratic ⇒ **every connected correlator of order `n≥3` vanishes identically** (Wick). Hence
```
χ^{(n)}_{connected} ≡ 0  (n≥3),   ⇒   𝔑_IR^{(n≥3)} ≡ 0,   for ALL m (gapped AND gapless).
```
In particular `𝔑^(3)(free massless) = 0 = ` the gapped value. The free massless field is **gapless** yet reads `𝔑^(3)=0`: a **false negative** — the detector labels the canonical gapless theory as gapped. This is independent of whether `𝔑` is a **frequency-IR** or **spatial-IR** (the `Λu` kernel) non-analyticity, because `χ^{(3)}` is identically zero before any IR analysis.

**Numerical confirmation `[V-here]` (Parts 1–2).**
- Part 1 (3D free field sampled by FFT diagonalisation, `N=2×10⁴`): as `m: 1.0→0.05`, the 2-pt variance grows (`0.170→1.001`, correlation length `ξ=1/m→∞`, i.e. **gapless in the IR**) while the connected 3-pt stays `|z|≤2.0` and the connected 4-pt cumulant `|z|≤1.3` — both consistent with **0**. The gaplessness lives **only** in the linear (KL) channel.
- Part 2 (tree-level frequency response `χ^{(3)}=g∏D(ωᵢ)`, `D=1/(ω²−m²+i0)`): gapped `m>0` is **smooth/analytic** (ratio `|χ³(10⁻³)|/|χ³(10⁻²)|=0.999`); gapless `m=0` **diverges** as `ω→0` (`|χ³|∝ω^{−6}`, fitted slope `−6.00`); **free `g=0` gives `χ^{(3)}≡0` for BOTH `m=1` and `m=0`** — the `g→0` limit is identical to the gapped value even though the spectrum stays gapless.

**Two structural corollaries `[C]`.**
- **`g`-continuity / fragility:** `𝔑^(3) ∝ g`. A *weakly* interacting gapless theory has `𝔑^(3)→0`, approaching the gapped value; the detector cannot separate "weakly-interacting gapless" from "gapped."
- **Operator / symmetry dependence `[V-here]` Part 5:** a **gapped** non-Gaussian theory has nonzero *equal-time* `χ^{(3)}_connected` *only if* the operator/symmetry allow it — a `Z₂`-symmetric `φ⁴` (gapped, interacting) gives `κ₃=0.0002±0.0008` (zero by symmetry), while explicit breaking (`h=0.4`) gives `κ₃=0.0136±0.0008` (16σ). The cubic correlator measures **interaction + order parameter + symmetry**, not the gap.

**Verdict on the nonlinear horn `[T,V-here]`.** `𝔑^(3)≠0` requires *(gapless)* **AND** *(non-Gaussian)* **AND** *(operator couples with non-analytic IR)*. The free-gapless case fails the middle condition outright. **The nonlinear discriminator is broken by free-field Gaussianity.**

---

## §2 · Solvable-model results

| model | gap? | Gaussian? | linear `𝔑` (analytic@0) | nonlinear `𝔑^(3)` (IR) | acid test |
|---|---|---|---|---|---|
| (i) free massive scalar | **yes** | yes | analytic `(𝔑=0)` `[T]` | `0` (`χ³≡0`) `[T]` | ok |
| (ii) **free massless scalar** | **NO** | yes | **non-analytic** (`M_1=∞`) `[T,V]` | **`0`** (`χ³≡0`) `[T,V]` | **FAILS** |
| (iii) interacting **gapped** (`φ⁴`, broken) | yes | no | analytic `[T]` | `0` (support>gap) `[C]`; `χ³_{eq}≠0` `[V]` | ok |
| (iv) interacting **gapless** (`φ³` lead. order; critical TFIM) | **NO** | no | non-analytic `[T]` | `≠0 ∝ g` `[C,V]` | ok **but ∝ g** |

- **(i)** free massive: linear sees the gap (`M_n` finite); `𝔑^(3)=0`; gapped ⇒ consistent. `[T/V]`
- **(ii)** free massless — **the decisive counterexample**: gapless, yet `𝔑^(3)=0` (= gapped value). The information exists **only** in the linear=KL channel. `[T,V]`
- **(iii)** `φ⁴` broken phase: gapped + non-Gaussian ⇒ nonzero *equal-time* 3-pt (`κ₃=0.0136±0.0008`), but its **frequency response is analytic** (`0∉` above-gap support) ⇒ `𝔑^(3)=0`. Consistent, but shows the equal-time cubic correlator is itself gap-orthogonal. `[C,V]`
- **(iv)** interacting gapless: `𝔑^(3)≠0` (tree-level `|χ³|∝ω^{−6}`), **but** `∝g` ⇒ `→0` in the free limit while remaining gapless. `[C,V]`
- **2D YM / CP^{N−1}** (the other suggested solvable models) are **gapped** (`[R]` standard: 2D YM via Migdal heat-kernel; CP^{N−1} dynamically gapped for `N≥2`); as gapped interacting theories they give `𝔑^(3)=0` (analytic) — **consistent, not discriminating**. They cannot serve as the gapless test case; the free-massless saddle (ii) does, and it fails.

**Decisive answer to the task's question.** The case **`𝔑=0` but gapless** exists and is the *simplest possible* theory: the **free massless scalar** (`𝔑^(3)=0` by Gaussianity). A gapless theory with a *decoupled* free massless sector (or one whose order/symmetry zeroes the cubic) is a second. **⇒ detector invalid.**

---

## §3 · Anomaly test — does the trace anomaly pin `𝔑`?

`[T]` Pure YM: `θ^μ_μ = (β(g)/2g) trF²`; the dilaton/anomaly low-energy theorem fixes the **convergent** moment of the `trF²` spectral density, `M_0 = ∫ρ(s)/s² ds` (the coefficient of the `ω⁰` Taylor term of `G_R`, §1.1).

`[V-here]` (Part 4) For the two measures:
```
gapped δ(s−0.7):   M_0 = 1/m² = 2.0408      (finite, positive)
gapless witness:   M_0 = ∫ s²e^{−s/Λ}/s² ds = Λ = 1.0000  (finite, positive)
```
**Both finite, both positive, same functional form.** The witness `dν=s²e^{−s/Λ}ds` is *positive, gapless (support down to `s=0`), and saturates the convergent sum rule*. **The anomaly moment cannot tell them apart** — it is **blind to the gap**. The gap is encoded precisely in the **non-convergent** higher moments (`M_1=∞` for the witness, §1.1), i.e. in the IR non-analyticity — which the anomaly does **not** fix.

**Conclusion `[V-here]`.**
- `𝔑` **defined by / pinned to the anomaly moment** ⇒ **non-discriminating** (returns the same value gapped vs gapless). **Horn (b).**
- `𝔑` **defined as the full-response IR non-analyticity** ⇒ discriminating in principle, but **computing it = deciding that all higher moments are finite = the gap = circular**. **Horn (a).**

So the anomaly either **pins `𝔑` into uselessness** or leaves it **circular** — it never yields a computable non-circular discriminator. The candidate **inherits the anomaly's blindness** in exactly the computable version.

---

## §4 · Citation check (charter §2)

Spot-checked the six invoked physics references (context, not load-bearing — but verified):

| ref (as cited) | status | actual content / accuracy |
|---|---|---|
| Nandkishore, Choi, Kim, **PRR 3, 013254 (2021)** | **REAL, accurate** `[R]` | "Spectroscopic fingerprints of **gapped** quantum spin liquids…" — nonlinear (2D coherent) spectroscopy diagnoses **gapped** spin liquids. **Directly complicates the candidate: nonlinear response is used for GAPPED systems.** |
| Li, Oshikawa, Wan, **PRX 11, 031035 (2021)** | **REAL, accurate** `[R]` | "Photon echo from lensing of fractional excitations in Tomonaga-Luttinger spin liquid" — `χ^{(3)}` in a **gapless** Luttinger liquid. |
| Fava, Gopalakrishnan, Vasseur, Essler, Parameswaran, **PRL 131, 256505 (2023)** | **REAL, accurate** `[R]` | "Divergent nonlinear response from quasiparticle interactions" — the divergent-in-time term is **"absent for free bosonic excitations."** **Directly confirms the central defect: nonlinear response ≈ interaction/Gaussianity detector, not a gap detector.** |
| **Sun et al., Nature 2025** | **UNVERIFIED — flag** `[R-flag]` | No *Nature* (main) nonlinear-response paper by "Sun et al." found. Closest matches: *npj Quantum Inf.* 2025 (quantum algorithm for nonlinear response); assorted *Nature Communications* 2025. **Could not confirm; treat as unreliable.** |
| Kaib, Möller, Valentí, **PRR 8, L012052 (2026)** | **REAL, accurate** `[R]` | "Nonlinear spectroscopy as a magnon breakdown diagnosis and its efficient simulation" — third-order 2DCS via Lanczos. |
| Yokoyama, Tada, **PRB 113, 174439 (2026)** | **REAL, accurate** `[R]` | "Diagnosing energy gap in quantum spin liquids via **polarization amplitude**" — a gap diagnostic via the *twist-operator/polarization amplitude*, **not** nonlinear response. Relevant only as another gap probe. |

**Net:** 5/6 verified and accurately described; **1 (Sun et al. Nature 2025) unverifiable — flagged.** Crucially, the *verified* ones are **evidence against the candidate's premise**: Nandkishore–Choi–Kim use nonlinear response for **gapped** spin liquids, and Fava et al. state the nonlinear signature is **absent for free bosons**. The literature confirms nonlinear response tracks **interaction / gapped-spin-liquid structure**, **not** gapped-vs-gapless.

---

## §5 · Verdict and the precise resolution

**(a) or (b)?** `𝔑_IR^(n)` is **(a) a reformulation** on the linear horn (it *is* the KL spectrum; circular) and **free-vacuous** on the nonlinear horn (Gaussianity zeroes it for the free massless field — a false negative that mislabels a gapless theory as gapped). If instead pinned by the trace anomaly, it is **(b) computable but non-discriminating** (blind convergent moment). It is **not (c)**.

**The trilemma (sharpest statement).** A usable "6th species" must be *(i)* discriminating, *(ii)* computable/non-circular, *(iii)* correct on the free-gapless limit. The candidate embeds the gap information in exactly the two channels that **cannot** satisfy all three:
- the **connected-response non-analyticity** — Gaussian-blind (`𝔑^(3)(free massless)=0`) and, when nonzero, `∝ interaction` and operator/symmetry-dependent;
- the **convergent anomaly moment** — gap-blind (`M_0` finite & positive for gapped *and* gapless witness).

There is no third channel in the proposal. Hence `𝔑_IR^(n)` does not distinguish gapped from gapless in a way that survives the free-massless acid test, the `g→0` continuity, the operator/decoupling arbitrariness, and the anomaly blindness.

**What a genuine 6th species would require `[O]`.** A functional that is provably nonzero **iff** `m₁>0` **by construction** (not by spectrum lookup), that is **manifestly nonzero at the free-massless point** *or* — more honestly — that **correctly returns the gapped value there while still discriminating everywhere it matters** (a property the connected nonlinear response fails, because it also returns the gapped value for *any* weakly-interacting-or-decoupled gapless theory). The linear=KL channel carries the correct information but is the gap restated. **The information the gap detector needs is provably absent from the connected-nonlinear channel (Gaussianity) and from the anomaly moment (sum-rule blindness).** Recovering it requires the full IR spectrum — i.e. the gap — regardless of `n`. This is the precise sense in which the candidate is a reformulation wearing a nonlinear costume, not a new species.

**One honest caveat (no overreach).** This tests the *specific* functional `𝔑_IR^(3)` as defined (connected nth-order causal response, IR-non-analytic part, anomaly-fixed convergent moments). It is a **no-go for this construction and its immediate variants** (any connected-`n` response, any convergent-moment functional), **not** a proof that no IR-sensitive discriminator of the gap can exist. Per charter §2 (`∀class ≠ ∀possible`), the door on a genuinely different 6th species stays open — but its gap information must come from a channel that is **not** Gaussian-blind and **not** convergent-moment-blind, and the present proposal is neither.

---

### Files
- Witness: `D:\Dev\math\research\unconditional\w7\ir_response_test.py` (ran `py -3.12`, ~9 s, single job).
- Numbers: `D:\Dev\math\research\unconditional\w7\ir_response_test_results.json`.

### Refs (external, spot-checked §4)
Nandkishore–Choi–Kim PRR 3, 013254 (2021); Li–Oshikawa–Wan PRX 11, 031035 (2021); Fava et al. PRL 131, 256505 (2023); Kaib–Möller–Valentí PRR 8, L012052 (2026); Yokoyama–Tada PRB 113, 174439 (2026). *Sun et al. Nature 2025 — unverified.* Standard: Källén 1952; Lehmann 1954; trace-anomaly/dilaton WZ (Adler–Bardeen 1969; Migdal–Shifman-type low-energy theorems).
