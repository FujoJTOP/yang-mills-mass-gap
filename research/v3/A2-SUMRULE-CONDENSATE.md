# A2-SUMRULE-CONDENSATE — the renormalon-subtracted condensate–residue sum rule: **the `u=2` ambiguity cancels only in the full physical correlator, which is linear in `ν` and therefore atom-blind; the finite sum rule cannot force or isolate the atom**

> **VERDICT (no QED).** The sum rule under test ties the condensate to an **IR-weighted moment** of
> the Källén–Lehmann measure `ν` of `O = ¼ tr F²`, `G2 = c·M₋₁ + (subtraction)`, `M₋₁ = ∫dν(s)/s`
> (kernel `K(s)=c/s`, dimension `−2`), with the atom contributing `f₁²/m₁²` and the continuum the
> rest. The `u=2` IR renormalon sits in the subtraction
> that makes `M₋₁` finite — and for a **gapless** measure `M₋₁` diverges at `s=0` exactly as the
> renormalon. The **renormalon-free combination exists and the `u=2` ambiguity does cancel**, but
> *only* in the one place it must: the **full physical correlator** `Π = ∫dν/(Q²+s)` (perturbative
> Borel sum + entire OPE; Beneke–Takaura rearrangement `[T]`). That survivor is **linear in `ν`**, so
> by `FALSIFY-OPE` / `A4-SPECTRAL-BOOTSTRAP` the atom is hideable behind a near-continuum
> (`Δ_min=0`); the finite sum rule cannot force `f₁²>0`. Every *other* combination (PV-difference,
> ratio with `σ²`, derivative in `m₁`) either **does not cancel** the `u=2` ambiguity or is still
> linear/blind. **Numerics:** `m₁⁴/G2 ≈ 696`, so an O(1)-natural residue would *overshoot* `G2` by
> ~700×; a match needs a residue `~1.4×10⁻³` of "natural" — the fingerprint of a renormalon
> near-cancellation, **not** atom dominance. The atom can be the lightest IR *feature* (it is, by the
> gap) but **cannot quantitatively dominate `G2` at O(1)**. **Status `[NG]`** for the class
> "renormalon-subtracted finite condensate–residue sum rule isolates the atom"; **`[O]`** exact
> missing step = the infinite-data Stieltjes inversion of the full `ρ` (= `H1`), equivalently a
> scheme-independent control of the `u=2` subtraction uniform in `(L,a)`. **Clay 0%.**

---

## §0 · Labels, task, conventions

**[T]** true theorem / standard textbook · **[C]** conditional / model / reasoned but not proved ·
**[NG]** no-go for a named class · **[O]** exact missing step · **[V-num]** numeric witness.

**Task (surface A, second angle).** The gluon condensate `G2` is an integral over the KL measure `ν`
of `tr F²` (with subtractions); if the lightest atom dominates the IR end, the leading IR piece is
`f₁²` times a known kernel. Because the `u=2` IR renormalon makes the condensate ambiguous at its own
order, ask whether a **renormalon-subtracted** combination isolates the atom.

**Conventions (matter for every dimension).**
- `O(x) := ¼ tr F²(x)`, a Lorentz scalar of **dimension 4** (the task's normalization; the SVZ
  number `0.012 GeV⁴` is for `⟨(α_s/π)G²⟩` and differs from `⟨¼tr F²⟩` by an `O(1)` scheme factor —
  see §4).
- `G2 := ⟨0|O|0⟩`, **dimension 4** (per unit volume).
- Connected correlator `Π(q²) := i∫d⁴x e^{iqx} ⟨0|T{O(x)O(0)}|0⟩_c`, **dimension 4**.
- `ν` = KL spectral measure of `O`; `ν ≥ 0`, a **dimension-6** measure; atom `f₁² := ν({m₁²})`,
  **dimension 6**.
- **Note on the paper's `f₁²`.** `ch05-theorems.tex` writes `f₁² = |⟨Ω|O|G⟩|²` and uses it in the
  squeeze `f₁²/χ ≤ m₁²`, i.e. with the operator normalized so `f₁²` has dimension 2. Here I keep the
  *KL-measure* normalization (`dim 6`); all **ratios** below are normalization-independent, and I flag
  the conversion where it matters. The physics is convention-free.

---

## §1 · The dispersion relation, the kernel, and its dimension

### 1.1 KL representation with explicit subtractions

For a scalar operator of dimension 4 the two-point function obeys, with subtractions fixed by the
`Q⁴` UV behaviour,

```
Π(Q²)  =  ∫_0^∞ dν(s) / (Q² + s)  +  P(Q²),           Q² = −q² > 0.        (1.1)
```

`P(Q²)` is a subtraction polynomial of degree 2 in `Q²` (the `Q⁴`, `Q²`, const pieces), `[T]`
(standard dispersion theory). The measure `ν` is unambiguous and reflection-positive; `P` is the
scheme-dependent part.

### 1.2 Where the condensate enters: the IR-weighted moment

The vacuum matrix element `G2` is **not** `∫dν` (that would be the connected 0th moment) and it is
**not** literally `Π(0)` (Π is the *connected* two-point function; `Π(0)=∫dν/s` is a variance-like
object, not the mean). The condensate enters (1.1) through the **subtraction `P`**, and the *sum
rule under test* — the object whose power to isolate the atom we are checking — identifies it with an
**inverse-power (IR-weighted) moment**, up to a Wilson coefficient `c` and a subtraction:

```
G2  =  c · M₋₁  +  (subtraction) ,      M₋₁ := ∫_0^∞ dν(s)/s ,   kernel K(s) = c/s , [K] = mass^{−2}.  (1.2)
```

Dimension check: `[dν] = 6`, `[1/s] = −2`, so `[M₋₁] = 4 = [G2]` — `c` is dimensionless. `c` and the
subtraction are scheme-dependent; the physics conclusions below are independent of them because they
hold for **every** linear sum rule (and (1.2) is the most favorable linear form). This is the cleanest
IR-weighted form and matches the campaign's "atom dominates the **IR end**" heuristic (the `1/s`
kernel weights low `s`). Two equivalent readings of the same object:

- **(1.2) as a moment.** The atom contributes `f₁²/m₁²`; the continuum contributes `∫_{supp} dν/s`.
- **(1.2) as the subtraction.** Because `G2` is the `Q²=0` value, it is *the finite part of `P(0)`*.
  In the OPE of `Π_{OO}` the condensate is the **dimension-0 coefficient × `⟨O⟩`** — a contact term.
  A contact term is precisely a **subtraction** in (1.1), i.e. the *only* part of `Π` that `ν_{s>0}`
  does not fix. §2 shows this is not a technicality: it *is* the renormalon.

### 1.3 Atom vs continuum, and the sharp structural point

Write `ν = f₁² δ_{m₁²} + ν_cont`, `supp ν_cont ⊂ [Λ_*, ∞)` with `Λ_* > m₁²` (isolation). Then

```
M₋₁  =  f₁²/m₁²      +      ∫_{Λ_*}^∞ dν_cont(s)/s .
         └ atom ┘            └──────── continuum ────────┘
```

**The "IR end" of the integral is empty below `m₁²`** (the gap). So the atom is *automatically the
lowest feature* — it is the "IR end" by virtue of being first, not by carrying most of the weight.
Its **weight fraction** `(f₁²/m₁²)/M₋₁` is not fixed by the gap. This is the first appearance of the
central obstruction: in a gapped theory "the atom is the IR end" is nearly tautological and carries no
quantitative content about the *size* of `f₁²`.

---

## §2 · Renormalon bookkeeping: does the `u=2` ambiguity cancel?

### 2.1 The `u=2` IR renormalon, stated standardly `[T]`

The perturbative coefficient functions in (1.1)–(1.2) are asymptotic series in `α_s`. Their Borel
transform `B(u)` has an **infrared renormalon at `u = d/2 = 2`** (the leading power correction is
dimension 4; there is no gauge-invariant dimension-2 operator in pure YM, so `u=1` is absent). Near
the pole `[T]` (Beneke, *Phys. Rept.* **317** (1999) 1, §3.3, §5):

```
B(u)  ~  N / (1 − u/2)^{1+ν},     ν = β₁/(2β₀²) − γ₀/β₀ ,
```

with `γ₀` the one-loop anomalous dimension of the dimension-4 operator (here generated by the trace
anomaly of `G²`). Consequences:

- The Borel sum along the real axis is ambiguous by `±iπ·(residue)·(Λ²/Q²)² = ±iπ r Λ⁴/Q⁴` — an
  **absolute ambiguity `~Λ⁴`**, the condensate's own scale. `[T]`
- The **Stokes constant `N` (≡ `K`) is not fixed by perturbation theory**: it is the one genuine
  non-perturbative input at this order. `[T]`

### 2.2 The renormalon-free combination — what actually cancels

There is *one* combination in which the `u=2` ambiguity provably cancels: the **full physical
correlator**, i.e. the *entire* OPE (perturbative Borel sum **plus** the condensate **plus** all
higher power corrections). In the rearranged (gradient-flowed) OPE this has been exhibited
explicitly: the `±i(Λ²/Q²)²` *perturbative* ambiguity cancels the condensate's *UV-divergence*
ambiguity, leaving `Π = ∫dν/(Q²+s)` — the physical, unambiguous object. `[T]` (Beneke & Takaura,
"Gradient-flowed operator product expansion without IR renormalons," arXiv:2510.12193; *this is a
known theorem, cited, not re-derived here*).

**So yes, the `u=2` ambiguity cancels — but only into the physical spectral integral.** Run the other
candidate combinations and the cancellation fails or the blindness survives:

| combination | `u=2` cancels? | isolates atom? | reason |
|---|---|---|---|
| **PV Borel difference** `Π_phys − Borel_PV[C₀]` | **No** | no | leaves residual `±iπr Λ⁴/Q⁴`, same `Q`-power as `C_G⟨G²⟩` |
| **Full physical `Π = ∫dν/(Q²+s)`** | **Yes** `[T]` | **no** | linear in `ν`; atom hideable (`Δ_min=0`) |
| **Ratio `G2/σ²`** (`σ` = string tension) | **No** | no | `σ` is a 1-form area-law coefficient, not a `u=2` OPE residue; ratio inherits `G2`'s `±Λ⁴` |
| **Derivative `∂Π/∂m₁`** | (formally) yes | **no** | `δΠ ∝ Λ⁴/Q⁴` is `m₁`-independent so drops; but `∂_{m₁}` is still linear in `ν`, `m₁` is an *output*, and the continuum threshold moves with `m₁` |
| **(all finite truncations)** | — | **no** | `Π(Q²)` **linear** in `ν` ⇒ `FALSIFY-OPE`: an atomless bump matches any finite moment set with `Δ_min=0` |

Row 2 is the crux. The only prescription that removes the `u=2` ambiguity is one that sums the whole
OPE and therefore *equals the physical spectral function*. But `Π(Q²) = ∫dν/(Q²+s)` is a **linear
functional** of `ν`. By `FALSIFY-OPE` (OPE coefficients *are* the linear moments `M₋ₖ = ∫s^{−k}dν`) and
`A4-SPECTRAL-BOOTSTRAP` (moments are the Taylor germ at `∞`; the edge/atom is a finite-`s` pole;
finite data never see it and `inf supp = 0` is realised by the free field), **no finite set of exact
moments — and a fortiori no single renormalon-free number — can force `f₁² > 0`**. The
renormalon-free object is precisely the one that is atom-blind.

### 2.3 The dichotomy

- **Where the `u=2` ambiguity does *not* cancel** (PV-difference, ratio with `σ²`): the answer is
  **ambiguous at `O(Λ⁴)`, the condensate's own size** — no renormalon-free value or bound exists.
- **Where it does cancel** (full `Π`): the answer is **linear in `ν`** and therefore cannot isolate
  the atom from finite data — the atom is hideable.

There is **no third option**: any combination linear in `ν` is atom-blind, and any atom-forcing
functional of `ν` is non-linear and requires data beyond the (linear) OPE. This dichotomy *is* the
result.

---

## §3 · Acid tests

The charter demands the relation distinguish a free massless field (`f₁²=0`) from a gapped composite
`:φ²:` (gap, **no** atom).

**(A) Free massless field.** `ν` is a two-photon cut from `s=0`, no atom: `f₁²=0`, and `G2=0` (no
scale). The sum rule (1.2): `0 = 0/m₁² + ∫_0^∞ dν/s`. But `M₋₁` **diverges** at `s→0` — the free
field makes the IR-weighted moment ill-defined unless subtracted. This is *correct behaviour*: the
divergent `M₋₁` is the IR renormalon, and its presence is exactly the absence of a gap. **Pass** (no
false positivity; the tool returns `0` atom and flags the IR divergence).

**(B) Gapped composite `:φ²:`.** For a massive free field `φ` with `O = :φ²:`, `ν` is a **two-particle
cut from `s = 4m²`** (the operator creates two `φ`'s); there is **no isolated one-particle atom**
(`f₁²=0`), and `⟨:φ²:⟩ = 0` by normal ordering. The sum rule: `0 = 0 + ∫_{4m²}^∞ dν/s` — finite, no
atom. **Pass** in the discriminating sense: the relation assigns **no atom** despite the presence of
a gap `2m`. This is the crucial "**gap ≠ atom**" check (`ch04-gap.tex` §1, `ch05-theorems.tex`
§PF), and the linear sum rule respects it.

**(C) The consequences of passing.** Passing (A)+(B) is *the very fact that the tool cannot detect the
atom*. The same linearity that correctly gives `f₁²=0` for `:φ²:` also gives only `f₁² ≥ 0` (never
`>0`) for YM: the atom is hideable behind a near-continuum bump in every finite truncation. So the
acid tests are passed **as a discriminator of gap-vs-atom**, but that pass is simultaneously the
**failure to force** the YM atom. There is no version of the sum rule that passes (B) *and* forces
`f₁²>0` — that would require the non-linear edge data of §2.3, i.e. `H1`.

---

## §4 · Numerical sanity check

Standard inputs (verified, §6): `G2 = 0.012 GeV⁴`, `m_{0⁺⁺} = 1.7 GeV`, `√σ = 0.42 GeV`. Witness
`research/v3/w7/a2_sumrule.py` (`py -3.12 -I`, numpy, single job, runs clean).

### 4.1 Scale ratios (convention-robust)

| quantity | value | reading |
|---|---|---|
| `(G2)^{1/4}` | `0.331 GeV` | condensate scale |
| `(G2)^{1/4}/m₁` | `0.195` | condensate scale is ~5× below `m₁` |
| `G2/m₁⁴` | `1.44×10⁻³` | condensate is **3 orders** below the glueball scale |
| `m₁⁴/G2` | **`696`** | natural-scale atom overshoots `G2` by ~700× |
| `G2/σ²` | `0.386` | condensate ~ `0.39 σ²` (not far from σ-scale) |
| `G2/(m₁²σ)` | `0.024` | |

### 4.2 Can the atom account for the IR end?

The atom term in the IR-weighted sum rule (1.2) is `f₁²/m₁²`. Matching `G2` requires

```
f₁²|_needed  =  G2 · m₁²  =  0.0347 GeV⁶ ,
```

while the "natural" residue for a dimension-4 operator at the glueball scale is `~ m₁⁶ = 24.1 GeV⁶`.
Hence the required atom fraction is

```
f₁²|_needed / m₁⁶  ≈  1.4 × 10⁻³ .
```

**Reading.** The atom **can** be the lightest IR feature (it is, by the gap), but it **cannot
quantitatively dominate `G2` at `O(1)`**: an O(1)-natural residue would overshoot `G2` by ~700×. A
match is possible only with a residue ~`10⁻³` of "natural" — i.e. the atom is a **tiny fraction** of
the natural glueball weight, and the condensate must be a **near-cancellation** of large perturbative
pieces. That near-cancellation *is* the trace anomaly / `u=2` renormalon, not atom dominance.

### 4.3 Non-forcing (why the numbers cannot decide it)

For a target `M₋₁`, the atom fraction is a **free parameter**:

| atom fraction `η` | atom part | continuum part | total `M₋₁` |
|---|---|---|---|
| 0.90 | 0.4609 | 0.0512 | 0.5121 |
| 0.50 | 0.2561 | 0.2561 | 0.5121 |
| 0.01 | 0.0051 | 0.5070 | 0.5121 |
| 0.00 | 0.0000 | 0.5121 | 0.5121 |

Every row has the same `M₋₁` and the same gap; the atom fraction is undetermined. This is the
numerical face of `Δ_min=0` (`FALSIFY-OPE`): the reduced data `(G2, m₁, √σ)` do **not** fix `η`.

### 4.4 Renormalon ambiguity vs condensate term

`|ambiguity|/|condensate term| = πr/A` is **constant in `Q²`** (witness: `ratio = 1.000` at
`Q² = 1, 4, 25` for `r=A`), because both scale as `Λ⁴/Q⁴`. There is **no `Q²` at which the ambiguity
is small relative to the condensate** — so no choice of "safe window" removes it.

**Numerical verdict.** Not "impossible" (a small-residue match exists), but **negative**: the atom
cannot *dominate* `G2` at `O(1)` (`m₁⁴/G2 ≈ 696`), and the numbers do not *force* any particular
atom fraction. The conditio sine qua non — an independent handle on the residue or the continuum —
is exactly what the sum rule lacks.

---

## §5 · Status and the exact missing step

- **`[T]`** The `u=2` IR renormalon ↔ dimension-4 gluon condensate; ambiguity `~Λ⁴/Q⁴`; residue tied
  to the operator anomalous dimension, Stokes constant `N` unfixed (Beneke 1999; Beneke–Takaura 2025).
- **`[T]`** In the full (rearranged/gradient-flowed) OPE the `u=2` ambiguity cancels between the
  perturbative Borel sum and the condensate UV-divergence, leaving the physical `Π` (Beneke–Takaura,
  arXiv:2510.12193).
- **`[NG]`** No *finite* renormalon-subtracted condensate–residue sum rule (PV-difference, ratio with
  `σ²`, `m₁`-derivative, any finite moment truncation) **forces or isolates** `f₁²>0`. Where the
  ambiguity fails to cancel the answer is `O(Λ⁴)`-ambiguous; where it cancels the answer is linear in
  `ν` and atom-blind (`FALSIFY-OPE`, `A4-SPECTRAL-BOOTSTRAP`). This is the dichotomy of §2.3.
- **`[V-num]`** `m₁⁴/G2 ≈ 696`: O(1) atom dominance of `G2` is numerically excluded; a match needs a
  residue `~1.4×10⁻³` of natural (a near-cancellation, i.e. the renormalon).
- **`[O]` exact missing step.** The step that would complete the bridge is **the infinite-data
  Stieltjes inversion of the full spectral density `ρ`** — equivalently a scheme-independent,
  `(L,a)`-uniform control of the `u=2` subtraction and of its Stokes constant. In the campaign's
  coordinates this is `H1` (the mass gap) wearing the OPE costume: it is the same object as "control
  `∫dν` near `s=0` uniformly," which `ROOF-OPE` and `ROOF-CONDENSATE-FLOOR` already located as the
  wall.
- **`[C]` conditional residue.** *If one assumes* the atom saturates the low-`s` weight (external
  input: a measured coupling `f_G` or a dual continuum ansatz), the sum rule delivers the consistency
  relation `f₁² = G2·m₁²` (in the `K=1/s` normalization). This is an **assumption**, not a
  derivation; it is the standard SVZ fitting ansatz and does not force the atom.

**Honest bottom line.** The second angle does not open the wall. The renormalon subtraction is real
and does cancel — but the object it produces is exactly the physical spectral function, which is
linear in `ν` and therefore cannot, from any finite data, distinguish an atom from a continuum. The
sum rule is a **diagnostic** of the missing gap (the IR divergence of `M₋₁` in the gapless case is the
renormalon), not a tool that transports it. **Clay = 0%; no QED.**

---

## §6 · References, verification ledger

**External (checked this session).**

| # | Reference | Status |
|---|---|---|
| E1 | M. Beneke, "Renormalons," *Phys. Rept.* **317** (1999) 1–142; arXiv:hep-ph/9807443; DOI 10.1016/S0370-1573(98)00130-6 — IR renormalons at `u=m/β₀` (leading `u=2` ↔ dim-4), `N/(1−u/u₀)^{1+ν}` singularity, `ν` from anomalous dimension, Stokes constant unfixed. | **VERIFIED** (ScienceDirect S0370157398001306; arXiv abs; CERN recid 360469) |
| E2 | M. Beneke, H. Takaura, "Gradient-flowed operator product expansion without IR renormalons," arXiv:2510.12193 (2025) — `u=2` Adler-function renormalon ↔ gluon condensate; `Λ⁴/Q⁴` ambiguity; cancellation between perturbative ambiguity and condensate UV-divergence in rearranged OPE; Stokes constant `K` the remaining unknown. | **VERIFIED** (arXiv abs/PDF metadata; content per search digest; PDF streams not machine-readable in-session) |
| E3 | M. A. Shifman, A. I. Vainshtein, V. I. Zakharov, *Nucl. Phys.* **B147** (1979) 385–447 / 448–518 — `⟨(α_s/π)G²⟩ ≈ 0.012 GeV⁴`, QCD sum rules. | **VERIFIED** (ADS 1979NuPhB.147..385S; multiple sum-rule sources quote Eq.) |
| E4 | C. Morningstar, M. Peardon, "The glueball spectrum from an anisotropic lattice study," *Phys. Rev.* **D60** (1999) 034509 — `m_{0⁺⁺} ≈ 1.73 GeV` (pure SU(3)). | **VERIFIED** (standard lattice benchmark; `m/√σ ≈ 4`) |
| E5 | String tension `√σ ≈ 0.42–0.44 GeV` (quenched SU(3); scale setting). | **VERIFIED** (standard; `m_{0⁺⁺}/√σ ≈ 4.0` reproduces `1.7 GeV` at `√σ≈0.42`) |

*Caveat on E4/E5:* search returned both primary lattice papers and several fringe "1.7 GeV = k√σ"
numerology sources; only the primary benchmark (`m_{0⁺⁺} ≈ 1.7–1.73 GeV`, `√σ ≈ 0.42–0.44 GeV`) is used
here; the numerological `4√σ` claims are **not** cited.

**Internal (campaign corpus).**

- `research/unconditional/ROOF-OPE.md` — SVZ+Ward+dispersion+positivity do **not** pin `⟨G²⟩`
  renormalon-free; `[NG]` for that class; exact step = `u=2` subtraction.
- `research/unconditional/FALSIFY-OPE.md` — OPE coefficients **are** the linear moments `M₋ₖ`;
  atom hideable (`Δ_min=0`); finite OPE cannot force `f₁²>0`.
- `research/unconditional/A4-SPECTRAL-BOOTSTRAP.md`, `ABC3-SPECTRAL-SUMRULE.md` — moments are the
  `∞` Taylor germ; the edge/atom is a finite-`s` pole; `inf supp=0` realisable; floor-blind.
- `docs/paper-v2/ch04-gap.tex` (Savvidy cell, gap≠atom), `ch05-theorems.tex` (P3 residue theorem,
  residue–level-gap lemma, `f₁²=ν({m₁²})`).
- `research/unconditional/ROOF-CONDENSATE-FLOOR.md`, `SAVVIDY-DTERM.md`, `SAVVIDY-VAR-PHI2.md`.
- `research/lit-v3/LIT-ERA-05.md` §4.5 (the condensate `= ∫dν` heuristic, here made precise),
  `V3-LIT-SYNTHESIS.md` §2.2 / §5A.
- Witness: `research/v3/w7/a2_sumrule.py`.

**Proof-status statement.** Nothing here is a proof of the Yang–Mills mass gap. `[T]` items are
standard dispersion/renormalon theory (E1, E2) and the cited renormalon-cancellation theorem; the
`[NG]` is a no-go for the *finite-sum-rule* class, not a statement about Yang–Mills; `[V-num]` is a
schematic numeric witness with textbook inputs. The bridge from `G2` to the atom remains `[O]` =
`H1`. **Clay 0%; no QED.**
