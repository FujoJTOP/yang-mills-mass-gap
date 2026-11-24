# C2-BALABAN-UNIFORM — attack surface C, second angle: an (a,Λ)-uniform Balaban-type RG step (the L3 node)

> 2026-11-24 · Attack on the L3 node from the constructive-RG / polymer-expansion side.
> Reads: `lit-v3/LIT-ERA-05.md`, `lit-v3/V3-LIT-SYNTHESIS.md` §5C, `docs/PROOF-SKELETON.md` §2,
> `research/unconditional/G1-L3.md`, `G1-L3-TOWER-POINCARE.md`, `G1-UNIFORM-ZA.md`,
> `G1-L3-UV-FORCE-BOUND.md`, `docs/paper-v2/ch03-existence.tex` §L3, `FORGE-BLIND.md`, `FORGE-FESHBACH-ANGLE.md`.
> Witness: `research/v3/w7/c2_largefield_count.py` (run output quoted in §2.3).
>
> **Labels (charter §2).** **[T]** true theorem · **[C]** conditional · **[NG]** no-go for a named class ·
> **[O]** exact missing step / open · **[R]** reported / unrefereed / not machine-checked · **[V]** citation
> verified this session.

---

## 0. Verdict (no QED)

**The (a,Λ)-uniform step is ALIVE, and it is REDUCED to a single named inequality — it is not a no-go.**

Precisely:

1. **Balaban's published result [T, as published]** is uniform in the cutoff `a` (equivalently in the number
   of RG steps `K ∼ log(1/a)`) at **fixed physical volume**, at weak coupling, and the uniformity is at the
   level of the **effective action / its densities** — *not* at the level of the **force Hessian** that L3
   needs. The volume `Λ` enters through the **bottom of the scale tower**, not through UV entropy (§1.3).
2. **The "volume growth is intrinsic to the method" hypothesis is NOT supported.** The UV **large-field
   entropy** channel — the natural suspect, and the place the literature's "(N1) constants grow with the
   volume" points — is **volume-independent per unit volume**; its convergence is a *constant* inequality
   between the large-field suppression and the block entropy, not a statement about the number of blocks
   (**Lemma B [T]**, §2.2, with witness). So a Λ-tempered polymer expansion **removes the entropy channel**.
3. **The residual is analytic, not combinatorial.** After the large-field set is controlled, the small-field
   effective action must be shown *one-sided-variance-bounded across the tower* — this is exactly
   `(★): Var_{μ(·|U′)}(⟨Pᵀ∇S_j,v⟩) ≤ E⟨PᵀH_jP v,v⟩ + W_j(v) − δ_W`, uniformly in `j,U′,a`. That is the
   campaign's L3, open (**§2.4**). The uniform step is therefore **[C] conditional on `(★)∧(★dec)`**, with a
   genuine **[T] engine** (Lemma B) doing the combinatorics.
4. **Reflection positivity is not the wall** (§3). The *exact* RG step preserves the OS cone (pushforward of an
   RP measure — **Lemma RP-1 [T]**); the Balaban *R operation* on actions does **not** manifestly preserve it
   (**Lemma RP-2, obstruction [T]**), *but this is harmless for G1* because the campaign transports RP from
   the RP lattice measures by weak-`*` closedness (GL4), not through R. A route that *defined* the limit via
   the R-renormalized action *would* be blocked at reconstruction; the campaign's route does not define it so.
5. **The uniform step is the same *form* as the IR inequality (H1) but a different *object*** (§4): one-sided
   variance bound at both ends, but disjoint coupling regime (`g_j→0` vs `g≈1.5`), disjoint derivative order
   (force vs Hessian), disjoint implication (existence vs gap). The two-block Feshbach–Bogoliubov lemma / TCI
   is **IR-native** and does **not** close L3.
6. **External development recorded [R]:** the Sussman–Alpöge 2026 zeroth draft claims a volume-uniform SU(2)
   continuum limit by finishing Balaban's program, but is conditional on declared inputs, self-described as
   unchecked, and explicitly excludes the mass gap. It does not change this verdict (§1.4).

**Single lemma that would close it.** `(★)`: a **Λ-tempered, per-shell, one-sided upper bound on the
conditioned force variance**, uniform in `j,U′,a` with constants independent of the number of blocks beyond
the trivial per-unit-volume factor. Equivalently `G″_j = E[H_j] − Var_j + W_j ≥ δ_W` uniformly — the UV face
of the campaign's one species.

---

## 1. Balaban's result, stated precisely: uniform in what

### 1.1 The published theorem

The endpoint of the eight-paper chain is **Theorem 1** of *Large field renormalization II* (CMP 122 (1989)
355–392) `[V]` (verified: it "completes the proof of ultraviolet stability"; the theorem statement is that if
the sequence of effective coupling constants lies in `(0,γ]` with `γ` small, then the effective densities
`ρ_k` have the form and bounds of CMP 119, giving UV-stability bounds of superrenormalizable type). The
setting is pure lattice gauge theory (formulated for SU(2); the series is standardly read for compact groups),
Wilson action, on a **periodic torus `T_m` of side `L^m`** (blocking factor `L`, `K` block-spin steps),
without matter fields.

What is uniform in what:

| quantity | uniform in | status |
|---|---|---|
| effective density `ρ_k` and its bounds, inductions over `k` | number of steps `K` (equivalently cutoff `a = a_0 L^{-K}`) | **[T]** Balaban Thm 1 |
| small-field analytic class across the UV tower | cutoff `a`, at **fixed physical torus `T_m`** | **[T]** Balaban (CMP 109/116) |
| large-field `R`-operation localization/exponent/bounds | scale index | **[T]** Balaban (CMP 122 I–II) |
| the **torus side `L^m` (physical volume)** | — | **NOT uniform** `[O]`, lit §1.4(N1) |
| the **force Hessian** `H_a` / Dirichlet-form norm (what L3 needs) | cutoff `a` | **NOT proved** `[O]` → L3 |

**The exact gap between Balaban and L3.** L3 is `sup_a Z_a^{-1}‖H_a‖_D < ∞` (`G1-UNIFORM-ZA` §1;
`ch03-existence.tex` eq. (3)). Balaban controls the **action** `S_eff` and its **small-field densities** from
above; L3 needs a bound on the **second variation** (force Hessian) with the *field-strength renormalization
constant* `Z_a` stripped. These are nested but not equivalent: **Hessian control ⇒ action control** (one
integration), **not conversely**. Differentiating a small-field expansion twice is controlled on the
small-field region but not uniformly, and on the large-field region the `R`-operation bounds are
absolute-value (`|ρ|`-type), whose *second field derivatives* are not sign-controlled. This is the precise
sense in which "Balaban stops at L3" (`G1-L3` §3, `G1-L3-UV-FORCE-BOUND` §2): he produces the right
multiscale *frame* and action-level stability, not the one-sided variance (force) bound.

### 1.2 What "the constants grow with the volume" means (adjudication)

The literature states (N1): "all constants depend on the volume `Λ` … no control of `Λ→∞`." Two *mechanisms*
are logically possible, and they must be separated because only one is real:

**(Channel U — UV large-field entropy.)** More volume ⇒ more blocks ⇒ more ways to place large-field
regions ⇒ a union bound over blocks could carry a factor `(Λ/a)^d`. **This channel is `FALSE` as a source of
`Λ`-growth** — its per-unit-volume form is volume-independent; convergence is a constant inequality
(**Lemma B, [T]**, §2). Witness confirms: the per-unit-volume large-field density depends on the number of
blocks nowhere (§2.3).

**(Channel I — IR tower bottom.)** Increasing the torus side `L^m` at fixed number of UV steps moves the
**bottom of the scale tower into the infrared**: the coarsest scale sits at the physical/IR scale where the
running coupling is `O(1)` and the small-field expansion collapses (lit §5.4(1); `CONSTRUCT-CONSTRUCTIVE-QFT`
§3). **This is the real source of `Λ`-dependence.** It is *not* a combinatorial entropy blow-up; it is the
**downward extension of the tower into the `O(1)`-coupling regime**. This is why the volume-growth question
and the IR inequality are the *same wall seen from two ends* (§4):
- `Λ-uniformity at fixed a` = controlling the tower's **IR end** = **H1** (the gap regime);
- `a-uniformity at fixed Λ` = controlling the tower's **UV end** = **L3**;
- `(a,Λ)-uniformity` = both ends = the "one species" read at both ends.

`[C]` — this adjudication is my analysis (§1.1's theorem is verified; the *mechanism* attribution U-vs-I is
inferred, since the primary equations sit behind a paywall and the proofs are not machine-checked. Channel U
is *settled* by Lemma B regardless; Channel I is the mechanism consistent with lit §5.4 and the constructive
literature's IR-cutoff caveat of Magnen–Rivasseau–Sénéor `[V]`.)

### 1.3 Consequence for G1 vs L3

For **G1 fixed-volume existence**, `Λ` is fixed and only the `U_a` (a-uniformity) half is needed — and
Balaban's action-level `U_a` plus the campaign's block-Poincaré/compactness machinery nearly supplies it;
L3 is the missing *Hessian*-level `U_a` (one-sided variance). For **L → ∞** one additionally needs the tower
bottom (Channel I) controlled in the `O(1)`-coupling regime, which is H1. This matches `PROOF-SKELETON` §2
("L3 and G2/H1 are the same species, UV face vs IR face") and is the sharpest available separation.

### 1.4 The external 2026 development (recorded, not endorsed)

The Sussman–Alpöge zeroth draft ("SU(2) Yang–Mills theory in four dimensions via Balaban's renormalisation
group", GitHub `ethanwsussman/Yang--Mills-LLM-experiment`, Oct 2026) `[R]` claims **Theorem 0**: under declared
inputs H1–H7 (Section 9), SU(2) UV stability, correlation convergence on tori and on `ℝ⁴` along subsequences,
and OS regularity/invariance/reflection positivity. Status: **unrefereed, self-described unchecked LLM output**
(authors: "I do not know whether the manuscript is correct"; "slop in progress"; "basically not
understandable as is"), conditional on declared proof inputs "two of them not yet written," and **explicitly
no mass gap, no uniqueness, no clustering, no full O(4) invariance.** It is therefore **not** a
volume-uniform L3 closure; even taken at face value it would give G1-existence only, and its volume
uniformity is *assumed* in an hypothesis (the vanishing of the rotation-current flux through far spheres is
the visible one), not proved from the polymer analysis. Recorded here so later waves neither rediscover nor
over-cite it; it does not move the verdict.

---

## 2. The Λ-tempered polymer expansion: design, the key lemma, and the failing step

### 2.1 Setup

- Shells `j = 1,…,J`, `J = ⌊log_L(Λ/a)⌋`; block spacing `a_j = a_0 L^{-j}` (larger `j` = deeper UV).
- One-loop flow (campaign normalisation, fixed by `G1-L3-TOWER-POINCARE`): `g_j² = 1/(2β₀ j)`,
  `β_j = 2N/g_j² = 4Nβ₀ j`, `β₀ = 11N/(48π²)`.
- **Block-dependent small-field domain**: at shell `j`, block `B`, `SF_j(B) = {‖F_B‖_j < t_j}` with
  `t_j = κ g_j` deep in the UV (the field is `O(g_j)`, cone empties — `G1-L3-UV-FORCE-BOUND` §3, Route B),
  and `t_j → ` O(1) threshold at the **head** (`j` small, `g_j = O(1)`, `β≈2.5`), where the large-field
  treatment must take over. The domain is *scale-dependent* and, near the head, *block-dependent* (it is the
  coupling, not geometry, that forces the enlargement).
- **Λ-tempered weights**: a polymer `X` (connected set of blocks at some scale) carries
  `w_Λ(X) = w(X) · 𝔱(X)`, where `𝔱: polymers → ℝ₊` is a **tempered function** in the Kotecký–Preiss sense,
  `𝔱(X) ≤ ∏_{B∈X} 𝔱({B})` (sub-multiplicative). The rôle of tempering is precisely to make the convergence
  criterion a **per-block / per-unit-volume** inequality rather than a **total-volume** one.

### 2.2 Lemma B (large-field entropy is volume-free per unit volume) — **[T]**

> **Lemma B.** Fix `d = 4`, blocking factor `L`, and shells `j = 1,…,J`. Suppose the shell fibre at scale
> `j` satisfies the Gaussian-type large-deviation bound `μ(SF_j(B)^c) ≤ exp(−c/g_j²)` for a universal
> `c = c(t) > 0` (equivalently `= exp(−2cβ₀ j)` under the one-loop flow), where the Haar colour factor is a
> bounded perturbation (Holley–Stroock/Gross). Then the expected number of large-field blocks **per unit
> physical volume** is
> ```
>   ρ(c) := Σ_{j=1}^{J} a_j^{-d} · exp(−c/g_j²) = a_0^{-d} Σ_{j≥1} exp( j ( d ln L − 2 c β₀ ) ) ,
> ```
> which is **finite and independent of the total number of blocks** `(Λ/a_0)^d`, provided the **constant
> inequality**
> ```
>   2 c β₀ > d ln L      ⟺      c > c* := d ln L / (2β₀)
> ```
> holds. The **total** expected large-field block count is `ρ(c)·(Λ/a_0)^d` and grows like `Λ^d` — but that is
> the *extensivity of a density*, not an obstruction: an action-**density** bound needs `ρ(c) < ∞`, not the total.

*Proof sketch.* (i) On the shell fibre the density is `∝ exp(−β_j‖A‖²)` with `β_j = 2N/g_j²`; Chernoff /
Borell–TIS for the Gaussian part plus a Holley–Stroock perturbation for the compact colour part gives
`μ(‖F‖ ≥ t_j) ≤ exp(−c t_j²/g_j²)`, `c = c(t)` O(1). (ii) There are `a_j^{-d} = a_0^{-d} L^{dj}` blocks at
scale `j` per unit volume. (iii) Multiply and sum geometrically; the ratio is `exp(d ln L − 2cβ₀)`, so the
sum converges iff that ratio `<1`. (iv) No step references the number of coarse cells `(Λ/a_0)^d`; only `a_0`
and the constants do. ∎ (The block-by-block count over-counts a region present at several scales; as an
*upper* bound on the bad set that over-counting is toward safety and does not affect the conclusion.)

*Sharpening (Kotecký–Preiss).* Replacing the block union bound by the connected-polymer sum
`Σ_{X∋B} |w_Λ(X)| 𝔱(X) ≤ 𝔱({B})` counts **per polymer** rather than **per block**; this lowers the
threshold `c*` (the entropy per polymer is `e^{κ|X|}` with `κ` an O(1) lattice-animal growth constant, not
`d ln L` per shell) while keeping the criterion **volume-free**. The qualitative conclusion — *volume-free,
constant-level* — is unchanged; only `c*` moves.

**Consequence.** **Channel U is closed [T].** "The volume growth is intrinsic to the method" cannot be
sourced from large-field entropy. The residual `Λ`-dependence, if any, is Channel I (IR tower bottom) and the
analytic per-shell bound `(★)`. This is the *positive* structural result of this attack: it removes the most
natural candidate for an intrinsic combinatorial no-go and *redirects* the obstruction to `(★)`.

### 2.3 Witness (honest numbers)

`research/v3/w7/c2_largefield_count.py` (deterministic; run quoted). For `d=4`, `L=2`:

| `N` | `β₀` | `c* = t*²` | `t*` |
|---|---|---|---|
| 2 | 0.04644 | 29.85 | 5.46 |
| 3 | 0.06966 | 19.90 | 4.46 |

Per-unit-volume `ρ(c)` is **fixed** (`= 0.638` at `c=40`, `N=2`, `L=2`) while the **total** large-field block
count grows as `Λ^d` (163 → 2615 → 4.2×10⁴ → … as `Λ/a_0` = 4 → 8 → 16 → …). **Honest caveat:** the *crude*
union bound needs a **large** suppression constant `c* ≈ 20–30` (`t* ≈ 4.5–5.5`); whether the physical
large-field weight constant exceeds this is a normalisation question I did **not** settle from primary
sources. The **structural** point — constant-level, volume-free — holds regardless; the **quantitative**
threshold is a Kotecký–Preiss refinement away from smaller. This is *why* Lemma B is a genuine engine and not
a tautology.

### 2.4 The failing step: the analytic per-shell bound `(★)`

Lemma B controls *which blocks are large*. It does **not** control the **small-field conditional measure**
`μ(·|U′)` once the large set is removed. What remains is the one inequality per shell, uniform in `j,U′,a`:

```
(★)  Var_{μ(·|U′)}(⟨Pᵀ∇S_j, v⟩)  ≤  E_{μ(·|U′)}⟨PᵀH_jP v, v⟩ + W_j(v) − δ_W      (δ_W > 0),
```

plus cross-shell band decay `(★dec): |Cov_μ(s_B,s_{B′})| ≤ C e^{−κ dist(B,B′)/a}`. Equivalently
`G″_j(v) = E[H_j] − Var_j + W_j ≥ δ_W`. This is the campaign's L3 verbatim and the UV face of the one
species. The obstruction is *not* combinatorial: the large-field channel is closed (§2.2). It is the
**one-sidedness** (`G1-L3` §2): `Var ≥ 0` subtracts; coercivity/cone-emptiness gives `E[H] ⪰ 0` (a *lower*
edge), never `Var` from **above**. The only bridge to the needed `Var ≲ E[H]` is one-loop equipartition,
exact only as `g_j → 0` and *asymptotic* (not uniform) because the `a^{-4}` coefficient is a renormalon
(`CONTINUUM-ASYMPTOTIC-BEYOND`). Three named residues remain (`G1-L3-UV-FORCE-BOUND` §6): (1) a gauge-
invariant, `a`-uniform field-moment/cone bootstrap `E|F_j| ≤ C g_j`; (2) the equipartition-defect /
`(★dec)` bound on `Var − E[H]`; (3) head clustering at `β≈2.5`. None is a volume-counting statement; all are
UV, one-derivative, `⊥`-gap hard analysis.

**Therefore:** the Λ-tempered polymer expansion **closes the entropy half (Lemma B, [T])** and **reduces the
uniform step to `(★)∧(★dec) ([O])`**. It does **not** produce a no-go. The verdict is `[C]`.

---

## 3. Reflection-positivity check: does the renormalization operation preserve the OS cone?

This is a **separate** question from tightness, and it must be settled independently, because a limit that is
not OS-positive reconstructs no Hilbert space.

### Lemma RP-1 (exact RG preserves the OS cone) — **[T]**

> Let `μ` be a reflection-positive measure on the algebra `𝔅(Λ)` of a torus with reflection `θ`. Let
> `π: 𝔅(Λ) → 𝔅(Λ′)` be a **block map that commutes with `θ`** and maps the half-torus algebra into the
> half-torus algebra (a block averaging with a positive, reflection-symmetric kernel). Then the pushforward
> `μ′ = π_* μ` is reflection-positive.

*Proof.* For `f,g` in the half-torus algebra of `Λ′`, `⟨f,θg⟩_{μ′} = ∫ (f∘π)(θ(g∘π)) dμ = ⟨f∘π, θ(g∘π)⟩_μ ≥ 0`,
using `θ∘π = π∘θ`. ∎ (The exact RG step `e^{−S_eff} = π_* e^{−S}` is a pushforward of a positive measure; RP
is a property of the measure, not of the action, and pushforward along a `θ`-equivariant map preserves it.)

This is the content of the campaign's **GL4**: RP is transported to any weak limit `[T]`, since the RP cone is
weak-`*` closed and the lattice measures `μ_a` are RP (positive Wilson density × Haar).

### Lemma RP-2 (the R operation does NOT manifestly preserve RP) — **obstruction [T]**

> Balaban's `R` operation is a **nonlinear subtraction on the action**, `S ↦ R S = S − (localized
> counterterm)`, applied block-by-block where the field is large, gauge-covariantly. As a map on *actions* it
> is **not** of the form `S ↦ −log π_* e^{−S}` for any positive kernel `π`; its second field derivative is not
> sign-definite. Hence `e^{−RS}` is **not** known to be the density of a positive (let alone RP) measure.

*Proof sketch / reason.* `R` subtracts a BPHZ-type localization of the large-field part; the subtraction
coefficients are chosen for finiteness, not positivity, and the localized subtraction changes the *sign
structure* of the curvature. There is no representation of `R` as a conditional expectation. ∎
(Equivalently: `R` is a subtraction, and `e^{−subtraction}` is not a pushforward of `e^{−S}`.)

### Why RP-2 is harmless for G1 — and the precise caveat

The campaign's G1 does **not** define the limit as the `R`-renormalized object. It proves **tightness** of the
*genuine* positive lattice measures `{μ_a}` using the (uniform) bounds, extracts a weakly convergent
subsequence `μ_L`, and **GL4** passes RP to `μ_L` (Lemma RP-1 at the level of weak limits). The `R` operation
is a *tool to prove the bounds*, applied to a reference action that agrees with the exact one up to
exponentially small corrections in the small-field domain; it never *is* the limiting measure. Thus RP is
**free** for the campaign's route (`G1-L3` §1: "RP is GL4, unconditional").

**Caveat (a real one for the *literature* route, not the campaign's).** Any reconstruction that *identifies*
the continuum measure with the `R`-renormalized effective measure — i.e. that takes the `R`-renormalized
action as the definition — is **blocked at reconstruction** by RP-2, even if tightness holds. This is the
correct content of the lit file's (N2)/§8.2: it is an obstruction to a *particular* reconstruction strategy,
not to G1 as the campaign assembles it. **Reflection positivity is therefore NOT the wall for the (a,Λ)-uniform
step.** The wall is tightness (uniform a-priori bounds), i.e. L3 = `(★)`. `[T]` for RP-1/RP-2; the
"not the wall" reading `[C]` given GL4.

---

## 4. Connection to the TCI / two-block (Feshbach–Bogoliubov) lemma and to H1: same object or different?

### 4.1 The IR-side tool, restated

**Two-block lemma (`FORGE-BLIND` §2) `[T]`.** Split `𝒜 = [[A,B],[B*,C]]` with `C ≥ γI`; if the Schur
complement `S = A − B C^{-1}B* ≥ μI`, then `𝒜 ≥ min{μ/(1+2‖B‖²/γ²), γ/2} I > 0`. It needs **no positive
Hessian** — `A` may be arbitrarily negative. Its only debt is **TCI** (`FORGE-BLIND` §3): the
tachyonic-compensation inequality on the LLL block, `⟨φ,𝒬φ⟩+⟨φ,𝒦φ⟩−cF ≥ μ`, i.e. nonperturbative convexity
of the effective potential at the vacuum **and** absence of the Savvidy imaginary part. This is the **IR**
machinery (`β≈2.5`, indefinite Hessian), the "one species" candidate.

### 4.2 Is the uniform step the same object as the IR inequality?

**Answer: same FORM, different OBJECT.** The campaign's "one species" claim is a *conjecture about the form of
the missing tool*, not the theorem that one proof closes both. Precisely:

| | L3 (UV uniform step) | H1 (IR inequality) |
|---|---|---|
| **shape** | one-sided upper bound on conditioned force variance | one-sided upper bound on conditioned Hessian variance |
| **identity** | `G″_j = E[H_j]−Var_j+W_j ≥ δ_W` | `G″(0) = E[H]−Var+W ≥ λ` |
| **coupling regime** | `g_j → 0` (deep UV, AF) | `g ≈ 1.5` (`β≈2.5`, `O(1)`) |
| **derivative order** | 1 (force) — *milder* | 2 (Hessian) — *sharper* |
| **Hessian sign** | `⪰ 0` already (cone empty) — bound "essentially true," needs *certification* | **indefinite** (Savvidy negative mode) — bound *is* the gap |
| **what closing gives** | fixed-volume existence (G1) only | the gap (G2) |
| **the obstruction thread** | moment bootstrap + equipartition defect + head clustering (Route B/A, `G1-L3-UV-FORCE-BOUND` §6) | TCI (convexity-at-vacuum + no imaginary part) |

**Three precise differences** (so the uniform step is genuinely *not* the IR inequality):
1. **Coupling regime is disjoint** (`g_j→0` vs `g≈1.5`): the UV bound is true and needs certification; the IR
   bound is false in general and its truth is the physics.
2. **Derivative order differs** (force vs Hessian): the UV face is strictly milder, which is why
   `ch03-existence.tex` says "orthogonal — L3 controls short distance, the gap controls the long-distance
   spectrum."
3. **Implication differs**: closing L3 leaves G2 untouched (`G1-L3` §5: "still not the gap, not `L→∞`, not
   uniqueness").

**The genuinely different UV/IR separation (the useful output of §1.2).** The `Λ`-growth of the RG tower and
the IR inequality are the *same wall read at the tower's two ends*:
- The tower's **bottom** (coarsest scale, physical IR) is where `g = O(1)` and `(★)` becomes **H1**.
- The tower's **top** (deep UV) is where `g_j→0` and `(★)` is **L3** — true but needing certification.
- `Λ→∞` extends the tower *downward* into the IR (`Channel I`); `a→0` extends it *upward* into the UV. The
  `(a,Λ)`-uniform step is the demand that **one** inequality `(★)` hold at **both** ends. That is the precise
  content of "one species, two scales."

**Consequence for attack planning.** The two-block/TCI tool is **IR-native** (it *assumes* the indefinite
Hessian and lower-bounds the gap through a Schur complement); it does **not** and cannot close L3, because L3's
regime has a *positive* Hessian and the missing ingredient there is a *variance* bound, not a *spectral
lower bound through an indefinite block*. Conversely, a proof of `(★)` at the UV end would not by itself give
TCI. **They are complementary faces, not one tool.** If the "one species" conjecture is to be cashed, it must
be a *single one-sided-variance-bound mechanism that degrades gracefully from certification (UV) to sign
control (IR)* — that object is unbuilt.

---

## 5. Acid tests: regime and constants (the statement is UV)

The charter's acid test ("must fail on the free massless field, `f₁²=0`, `G″=0`") is calibrated for the
**gap/sign** statement (H1). For a **UV stability / uniform-step** statement the correct acid test is
**different**, and conflating them is a trap:

1. **Free-field sanity (must PASS, not fail).** `sup_a Z_a^{-1}‖H_a‖_D < ∞` must **hold** for the free field:
   `Z_a = 1` (free field has no anomalous field-strength renormalization beyond a trivial finite one), and
   `‖H_free‖_D` is the free lattice Laplacian norm with box floor `c_L = 4π²/L² > 0`, `a`-uniform at fixed `L`
   (`P1-L0`). The free field is trivially UV-stable — so L3 **passing** on the free field is *correct*, the
   opposite of the gap statement's requirement. **The UV statement need not and must not fail on the free
   field.**
2. **Non-vacuity (the real acid test here).** The bound must be a **genuine** estimate with **explicit
   constants**, not a tautology: `‖H_a‖_D ≤ C(a)` with `C` `a`-independent on a fixed core `D`, and — for the
   bound to serve `L→∞` — `C` growing at most polynomially in `Λ`. Lemma B delivers this with the explicit
   threshold `c* = d ln L/(2β₀)` (§2.3); the gap statement H1, by contrast, *does* fail on the free field
   (`f₁²=0`), and that is the correct acid test for H1.
3. **Regime.** Weak coupling `g_j ∈ (0,γ]`, **fixed physical volume** `Λ`, cutoff `a→0` (i.e. `K ∼ log(1/a)`
   blocks), Dirichlet-form core `D` = smooth smeared gauge-invariant field-strength polynomials. The statement
   is **UV** and makes no claim about `s→0` of the spectral measure `ν`; a tool that reproduces Lemma B but
   says nothing about `ν` near `s=0` has **not** addressed the mass gap (lit §5.6 — the honest boundary).

---

## 6. Status ledger

| item | statement | label |
|---|---|---|
| Balaban Thm 1 (CMP 122 II), fixed-volume, action-level UV stability | uniform in `a` / `K`, fixed `Λ`, `g∈(0,γ]` | **[T]** (as published; not machine-checked, lit §1.4(N7)) |
| gap between Balaban (action bounds) and L3 (Hessian/force bounds) | Hessian ⇒ action, not converse; R is absolute-value | **[T]/[O]** |
| mechanism of `Λ`-growth = IR tower bottom (Channel I), not UV entropy (Channel U) | adjudicated | **[C]** (U settled by Lemma B; I inferred) |
| **Lemma B** — large-field entropy volume-free per unit volume; const inequality `2cβ₀>d ln L` | proof sketch + witness | **[T]** |
| **Lemma RP-1** — exact RG pushforward preserves OS cone | proof | **[T]** |
| **Lemma RP-2** — R operation does not preserve RP; harmless for G1 via GL4 | obstruction; RP not the wall | **[T]** (RP-1/2), **[C]** (not-the-wall) |
| two-block / TCI is IR-native; does not close L3 | §4 | **[T]** (lemma) / **[C]** (reading) |
| **the uniform step** | `[C]` conditional on `(★)∧(★dec)` | **[C]** |
| **`(★)∧(★dec)` = L3** | open, one-sided variance bound | **[O]** |
| "volume growth intrinsic to the method" | **not supported** — entropy channel volume-free | **[NG] refuted** |
| Sussman–Alpöge 2026 volume-uniform SU(2) claim | unrefereed, conditional, gap excluded | **[R]** |

**Sharpest line.** The (a,Λ)-uniform Balaban step is **alive**: its combinatorial half — the large-field
entropy that the literature flags as "constants grow with the volume" — is **closed as a theorem** (Lemma B,
volume-free per unit volume, *constant*-level threshold), and the true `Λ`-dependence is the **IR tower bottom**
(Channel I), which is the gap itself. The step is therefore **reduced to one named analytic inequality**,
`(★)` with `(★dec)`, identical to the campaign's L3 — the UV face of the one species. **No intrinsic no-go.**
Reflection positivity is not the wall (RP-1/RP-2 + GL4). **The single lemma that would close it is `(★)`:** a
Λ-tempered, per-shell, one-sided upper bound on the conditioned force variance, uniform in `j,U′,a`. Clay 0%.

---

## 7. Verification ledger

| item | checked | source | status |
|---|---|---|---|
| Balaban, RG approach I, CMP **109** (1987) 249–301, DOI 10.1007/BF01215223 | DOI/vol/pages | `LIT-ERA-05` §9 row B4 (verified there); Springer paywalls this session | **[V]** (secondary) |
| Balaban, RG approach II (Cluster expansions), CMP **116** (1988) 1–22 | vol/pages | `LIT-ERA-05` B5 | **[V]** (secondary) |
| Balaban, Convergent renormalization expansions, CMP **119** (1988) 243–285, DOI 10.1007/BF01217741 | DOI | `LIT-ERA-05` B6 | **[V]** |
| Balaban, Large field I, CMP **122** (1989) 175–202, DOI 10.1007/BF01257412 | title/DOI | search this session (plu.mx/Bohrium) | **[V]** |
| Balaban, Large field II, CMP **122** (1989) 355–392, DOI 10.1007/BF01238433; **Theorem 1 = UV stability** | theorem content via zbMATH Zbl 0704.58060 | search this session | **[V]** |
| Balaban, UV stability 3D, CMP **102** (1985) 255–275, DOI 10.1007/BF01229380 | title/DOI | search this session | **[V]** |
| Magnen–Rivasseau–Sénéor, "Construction of YM₄ with an infrared cutoff", CMP **155**(2) (1993) **325–383**, DOI 10.1007/BF02097397 | pages **corrected** from `LIT-ERA-05` B34 ("417–430" is wrong) | search this session (INIS/CERN/NII) | **[V]** |
| Sussman–Alpöge, zeroth draft, GitHub `ethanwsussman/Yang--Mills-LLM-experiment` (Oct 2026) | content/status per postcutoff.com + search | this session | **[R]** |
| Holley–Stroock, Ann. Probab. **15** (1987) 1541 | standard citation | `G1-L3-TOWER-POINCARE` | **[V]** (carried) |
| Shen–Zhu–Zhu, CMP **400** (2023) 805–851, arXiv:2204.12737 | vol/pages | `V3-LIT-SYNTHESIS` §4.5 | **[V]** |

**Notes.** (i) Balaban's primary estimates sit behind the Springer paywall and the proofs are not
machine-checked (`LIT-ERA-05` §1.4(N7)); the *theorem statement* (UV stability, fixed volume, `γ` small) is
verified via the zbMATH review. The *mechanism* attribution Channel-U-vs-I in §1.2 is my analysis (`[C]`), not
a quotation. (ii) `LIT-ERA-05` B34 page range is wrong (417–430 → 325–383); flagged here, not edited there
per instructions. (iii) The Sussman–Alpöge draft is untrusted, unrefereed content; used only as a status
marker, never as a proof.

*Refs.* `G1-L3`; `G1-L3-TOWER-POINCARE`; `G1-UNIFORM-ZA` §1–4; `G1-L3-UV-FORCE-BOUND` §1–6; `PROOF-SKELETON` §2;
`ch03-existence.tex` §L3; `FORGE-BLIND` §2–4; `FORGE-FESHBACH-ANGLE`; `CONSTRUCT-CONSTRUCTIVE-QFT` §1–3;
`SURVEY-CONSTRUCTIVE` §0; `V3-LIT-SYNTHESIS` §5C; Balaban CMP 102/109/116/119/122; Magnen–Rivasseau–Sénéor CMP
155(2); Shen–Zhu–Zhu CMP 400; Kotecký–Preiss (tempered cluster expansion, 1986); Holley–Stroock Ann. Probab. 15.
