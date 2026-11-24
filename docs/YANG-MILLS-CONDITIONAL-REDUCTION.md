# A Conditional Reduction of the 4D Pure Yang–Mills Mass Gap

### and a record of what was actually proven

> **Status of this document.** This is a **conditional reduction and an audit**, not a proof of the Clay problem. The Clay mass gap for 4D pure SU(N) Yang–Mills remains open. The only defensible headline is: *the problem has been rewritten onto one precisely named, acid-tested object, and every honest route to it terminates there.*
>
> **Label legend used throughout.** Each substantive claim carries one of:
> **[T] theorem** — proved (any imported standard result is named); **[C] conditional theorem** — a proved implication from explicitly named hypotheses; **[N] no-go** — a theorem that a whole class of approaches cannot work; **[M] exact missing step** — a precisely stated open inequality, with no slack hidden; **[E] numerical evidence** — measured, not proved.
>
> **Provenance rule.** Every theorem below names the file that contains its proof. No claim is asserted without a source. Citations in the reference list appear verbatim in the corpus; where the corpus itself flags a bibliographic uncertainty, this is reproduced rather than silently repaired.

---

## Abstract

We report the outcome of a sustained campaign — Lean 4 formalization plus lattice numerics — aimed at the 4D pure SU(N) Yang–Mills mass gap. **We do not prove the mass gap.** We prove, instead, a reduction and a ledger.

The reduction has three parts. *(i) A rewriting of the problem.* The gap is equivalent to the statement that the Källén–Lehmann spectral measure `ν` of the gauge-invariant field-strength operator `tr F²` has its lightest atom `m₁²` at positive mass with positive weight `f₁² = ν({m₁²}) = |⟨0|tr F²|G⟩|²`, and that this atom is isolated. *(ii) An exact spectral squeeze.* For any such `ν` with susceptibility `χ = Var(a) = ∫ s⁻¹ dν`, one has `f₁²/χ ≤ inf supp ν ≤ M₀/χ`. It follows that a lower bound on the mass requires **both** a residue floor `f₁² ≥ r₀` and a susceptibility ceiling `χ ≤ 1/λ`; the often-quoted condition H1 (`χ ≤ 1/λ`) is **necessary but not sufficient**, because it is blind to a thin atom. The correct sufficient condition is the pair `(f₁² ≥ r₀) ∧ (χ ≤ 1/λ)`, which gives `m₁² ≥ r₀ λ`. *(iii) A convergence.* Every route we could construct — and we constructed on the order of ninety — ends at **one species**: a nonperturbative, non-circular, sign-controlling **upper bound on a conditioned variance**, appearing once in the ultraviolet (the L3 node of existence theory) and once in the infrared (H1). The same inequality, drawn at two scales, is the remaining gap between the campaign and Clay.

The ledger has a positive half and a negative half. On the positive side we record the genuine theorems: an unconditional strong-coupling **residue** theorem (`f₁² ≥ r₀(β) > 0` for `β < 1/48`, volume- and L-uniform — a statement about the atom weight I1, **not** its isolation I2); a kinematic **no-dissolution** theorem for the lightest `0⁺⁺`; the **non-negativity of the Haar geometric curvature** `W_geo ≥ 0` by Prékopa; a **gauge-invariant Perron–Frobenius** theorem on the physical Hilbert space; the **strong-coupling coarse Poincaré** lemma (Lemma SC); the abstract **Feshbach–Bogoliubov two-block lemma**; and the **exact variance decomposition** `G″ = S − D`. Six new tools pass the free-field acid test (they must fail on the massless free field, where `f₁² = 0` and `G″ = 0`). On the negative side we record the no-gos: the Schur-convex master obstruction, Goldstone uniqueness at zero mass only, the periodic, nested-cone, rigidity, Feshbach-angle and Krein obstructions.

Numerically, the residue is seen directly where it is resolvable: `f₁² > 0` at `β = 2.2` at seventeen sigma on SU(2) (`L = 6`); at `β = 1.0` it dissolves into the cutoff, exactly as the analytic strong-coupling picture predicts, not as a failure. A `k → 0` curvature probe on `L = 8` finds the coarse-mode `G″` positive, `W_geo`-carried, and volume-stable — but it is measured on a **gauge-dependent gluon mode**, not the `0⁺⁺` object, and its rigorous-style one-sided interval still includes zero.

The verdict is sober: **Clay is 0% — open.** The wall did not fall. It became transparent, ruled, and instrumented. The tool that would cross it — a nonperturbative machine transporting the strong-coupling gap and residue across the Savvidy crossing and taking `a → 0` — **does not exist**, and naming it precisely is the campaign's real result.

---

## 1. Introduction

### 1.1 The Clay statement, exactly

The Clay problem, as posed by Jaffe and Witten, is:

> Prove that for any compact simple gauge group `G`, a non-trivial quantum Yang–Mills theory exists on `ℝ⁴` and has a mass gap `Δ > 0`. Existence includes establishing axiomatic properties at least as strong as those cited in Streater–Wightman (1964) and Osterwalder–Schrader (1973, 1975). [Jaffe–Witten, *Quantum Yang–Mills Theory*, Clay Mathematics Institute]

Two features of this statement are habitually under-weighted and drive everything that follows. First, it demands a **uniform** gap: not merely that each finite lattice volume has a positive spectral gap (trivially true by compactness), but that the gap survives `L → ∞` and `a → 0` with a volume- and cutoff-uniform lower bound, and that the resulting theory clusters exponentially at a rate `m > 0`. Second, it demands that the gap be a feature of the **gauge-invariant, physical** theory — not of a gauge-fixed chart, and not of a perturbative expansion around a chosen background.

The campaign's operational form of the Clay conclusion, assembled and partially machine-verified in `docs/conditional-theorem.md`, is:

```
∃ μ on ℝ⁴, gauge-invariant, OS0–OS4;
   non-Gaussian, g_ren/Λ > 0;
   Δ/Λ ≥ min(δ_W, δ_S) > 0;
   |S(f τ_a g) − S(f) S(g)| ≤ C e^{−m‖a‖},  m ≥ c Δ > 0, all gauge-invariant local f, g;
   unique vacuum; area law σ = c_σ Λ² > 0; all constants group-uniform in (h^∨, dim G, dim f).
```

The transfer-matrix form of the same statement is the one used throughout this paper: for the lattice transfer matrix `T_a = e^{−H_a}` on the gauge-invariant Hilbert space `ℋ_a`, with vacuum `Ω`,

```
∃ δ > 0 :  inf spec( H_a |_{Ω^⊥} ) ≥ δ Λ    for all sufficiently small a.
```

### 1.2 The lock-core rewriting

The campaign works with a sharper, channel-specific rewriting of the gap which, while not literally equivalent to the Clay statement, is the correct object to attack. Let `ν` be the Källén–Lehmann spectral measure of the gauge-invariant local operator `O = tr F²` (more precisely `O_S = ¼ tr F_{μν}F^{μν}`, the scalar density; the normalization is immaterial). Reflection positivity gives `ν ≥ 0` and `supp ν ⊂ [m₁², ∞)` where `m₁` is the lightest mass in the `0⁺⁺` channel. The **atom** is

```
f₁² := ν({m₁²}) = |⟨0| tr F² |G⟩|² ,
```

the squared overlap of the vacuum with the lightest `0⁺⁺` glueball `G`. The **lock core** is: `f₁² > 0` and the atom is isolated (`Δ > 0` above it). The campaign summarises this as

```
Clay  ⟺  ¬(iii),   (iii) = "the YM measure is a gapless interacting 4D CFT".
```

This rewriting is useful because it names the exact analytic object (an atom weight in a positive measure at the bottom of a channel) and the exact analytic property it has: `f₁²` is **non-diagonal** (it is an off-diagonal matrix element squared), **upper-semicontinuous** in the coupling (`ν({x})` is usc), and requires a **second scale** `E₂` to certify isolation. Each of these properties will matter: usc means continuity can only make an atom *fall*, never rise; the need for `E₂` is the reason every bottom-functional tool is blind.

**Honesty checkpoint.** "Clay `⟺` atom" is *not* a theorem of equivalence. The campaign's own adversarial audit (`research/unconditional/REDUCTION-AUDIT.md`) establishes that the correct logical shape is a **diamond**, not a chain: a gapped theory need not have an atom at all (the massive free scalar with `O = :φ²:` has a gap and `f₁² = 0`), so `gap ⇒ atom` is false; the true direction is `isolated atom at inf supp ν ⇒ exponential clustering in that channel`, conditional on that channel carrying the lightest state. This paper therefore treats "the atom `f₁² > 0`" as the **target object**, not as a synonym for Clay. The gap between the two is itself a named missing step — the ordering input that the lightest state lies in `0⁺⁺` (the corpus's `IN-11`, a physical input, not derivable from symmetry).

### 1.3 The two halves

The campaign splits Clay into two halves and a channel:

- **G1 (existence).** 4D continuum existence with OS axioms. The fixed-volume construction is essentially complete; the one genuinely open UV input is a **uniform renormalization constant** `Z_a`, equivalently a cutoff-uniform bound on the force Hessian (the node called L3). L3 is purely ultraviolet and *orthogonal* to the gap — proving it would make fixed-volume existence unconditional but would not help G2.
- **G2 (gap).** The proof or transport of `f₁² > 0` (equivalently of the coarse curvature `G″(0) > 0` at `k → 0`). This is where the problem lives.
- **The physical channel.** Identifying the lightest state as `0⁺⁺` and showing it does not dissolve (the no-dissolution theorem). This half is, conditionally, complete.

The single deepest structural finding of the campaign is that **G1's L3 and G2's H1 are one species**: both are a one-sided *upper bound on a variance*. L3 is the UV force-variance bound; H1 is the IR conditioned-variance bound `Var_μ(⟨v,∇S⟩|0) ≤ W_geo`. The shapes are identical because `Var ≥ 0` can only supply a lower edge, while in the exact identity `G″ = E[H] − Var + W` the variance enters with a minus sign. The two faces differ only in scale (UV versus IR). This unification is stated in the corpus as a **structural identification and a conjecture**, not a theorem (`G2-ASSEMBLE-100.md` §7: "a structural identification … one closure would likely close both faces — but that is a conjecture, not a theorem"). This paper preserves that label.

### 1.4 What this paper is, and is not

This paper is **not** a proof of Clay, and it contains no step labelled QED. It is:

- a precise statement of the reduction, with every hypothesis displayed;
- a complete inventory of what was actually proved, labelled by type;
- a statement of the no-gos, each with the precise reason its class fails;
- the immune map: roughly ninety routes, each with the exact point at which it dies;
- a specification of the missing tool: five conditions, verbatim;
- an honest verdict.

The reader who wants only the bottom line should read §9. The reader who wants to check the work should follow the file names in §4.

### 1.5 Roadmap

§2 states the two halves and the node board with honest percentages. §3 gives the exact spectral squeeze and explains why H1 is necessary but not sufficient. §4 states the proven theorems with their hypotheses. §5 presents the six tools and each one's free-field acid test. §6 records the no-gos. §7 states the one-species reduction and the immune map. §8 specifies what the missing tool must satisfy. §9 gives the verdict. §10 collects verified references.

---

## 2. The two halves: G1 (existence) and G2 (gap), and the node board

### 2.1 G1 — existence

G1 asks for a 4D continuum SU(N) gauge theory satisfying OS0–OS4. The campaign's decomposition (`docs/PROOF-SKELETON.md` §2) is:

| Node | Object | Local % | Label |
|---|---|---|---|
| G1.0 | Fixed-volume assembly: L1 definition, L2 theorem, single-block `G″ ≥ δ_W`, GL3 compactness, GL4 reflection-positivity propagation | ~95% (all but L3) | [T] modulo L3 |
| G1.1 | Perturbative/combinatorial kernel: BRST dimension-4 classes reduce to `{Tr F², Tr F F̃}`; asymptotic freedom absorbs the marginal log; `Σ g_j⁴ = ζ(4)/(16β₀⁴) < ∞` | 100% | [T] closed |
| G1.2 | **Uniform `Z_a` = L3**: `sup_a Z_a⁻¹ ‖H_a‖_D < ∞` = L2 + block Poincaré + `Λ_UV` tower uniform | ~40% | [M] open |
| G1.3 | `a → 0` continuum via OS reconstruction, given L3 | ~70% | [C] on L3 |
| G1.4 | `L → ∞`, OS2 (continuous O(4)), unique vacuum OS4 (clustering) | 0% | [M] open, entangled with G2 |

The perturbative kernel G1.1 is a genuine closure: the marginal logarithm is not the obstruction (`Σ g_j⁴` converges), and the only surviving dimension-4 operators are `Tr F²` and `Tr F F̃`. G1.2 is the sole open UV node, and it is exactly a **single-sided force-variance bound** — the same species as G2's H1. G1.3 is conditional on G1.2 (OS0/OS1/OS3 are otherwise assembled; OS0 requires L3, OS1 and OS3 are proved). G1.4 (OS2, infinite volume, uniqueness) is open and, crucially, **OS4 is downstream of the gap**: a uniform gap gives clustering for free.

The campaign's honest aggregate: **G1 is 100% assembled modulo `{L3, G2}`**, with 18 input rows each classified as a theorem, a reduction to L3, the gap, or the same species — and no genuinely independent third open item. The master document records a coarser figure of "G1 ≈ 70%" (`docs/YANG-MILLS-MASTER.md` §9); that number is a summary of "perturbative kernel closed, L3 open," not a separate claim. We report both and do not average them.

### 2.2 G2 — the gap

G2 asks for `f₁² > 0` with isolation, equivalently `G″(0) ≥ λ` at `k → 0`, uniform in `(L, a)`. The decomposition (`research/unconditional/G2-ASSEMBLE-100.md`) is:

| Node | Object | Local % | Label |
|---|---|---|---|
| G2.0 | Identity layer: `G″(k→0) = E[H|0] − Var(∇S|0) = V_eff″ = m_eff²(k→0)`; Ward `G″·Var(a) = 1`; Schwinger–Dyson `E[H] = E[f²]` | 100% | [T] exact |
| G2.1 | Haar geometric term `W_geo = E[Φ″] − Var(Φ′) ≥ 0` (Prékopa), every mode including `k → 0`, every L | 100% | [T] |
| G2.2 | Shallow-field band `{2g‖F‖ ≤ N/12}` closed by `Φ″ ⪰ N/12` + Brascamp–Lieb | 100% | [T] |
| G2.3 | High-field cone (WSTC): cone variance mass `≤ (cβκ_max)² e^{−cβ(κ²−κ_typ²)N_□} → 0`, L-uniform | 100% | [T] |
| G2.4 | AF-UV deep cone: instability `≤ μ² e^{−cμ⁴/g⁴} → 0` | 100% | [T] |
| **G2.5** | **Open cell = typical field × high momentum**: `Var_μ(⟨v,∇S⟩|0) ≤ W_geo` | **0%** | **[M] = H1** |
| **G2.6** | **Named species**: negative association / anti-FKG `Cov(G, Φ′|0) ≤ 0` | **0%** | **[M] proven absent from nine classes** |

The region coverage of the `(field × momentum)` plane is about **90%**, closed by the theorems G2.1–G2.4. The remaining one cell — typical field, high momentum — is a region where all four handles fail *simultaneously*: the WSTC exponent vanishes when the field is typical, the AF-UV suppression has not begun at `g ≈ 1.5`, and the `N/12` floor is roughly ten times too small. That cell is **exactly** `Var(a_{k→0}) ≤ 1/λ`, with zero slack. The sign — the actual proof — is **0%**.

The G2 assembly is complete in the following precise sense: of 26 input rows, every row is either (a) a theorem or (b) a reduction to the single species; **there is no row (c)**, a genuinely independent open item. The curvature face H1 and the spectral-edge face H* (`f₁² ≥ r₀` with isolation I2) are one object in two vocabularies, tied by the residue–level-gap lemma (`I1 ⟺ I2`) and Ward (`G″ · Var(a) = 1`). The unification `H* = H1` is the conjecture; the dictionaries are theorems.

### 2.3 The dimensional correction

One correction must be recorded because it is easy to get wrong and it changes what "`G″ ≥ λ`" means (`docs/PROOF-SKELETON.md`, G1-CONTINUUM-SCALING). For the dimensionless Wilson action, the Hessian satisfies

```
G″_a = a² · m_phys² ,      m_lat = √G″_a = a · m_phys ,      Δ = lim m_lat / a .
```

Therefore a claim of the form "at fixed `a`, `G″_a ≥ λ` with `λ` a positive *constant*" is dimensionally wrong for the gap: it would give `m_phys = lim √λ / a = +∞`, the free-box UV scale, a lattice artifact rather than a particle. The correct target is a **ratio** `G″_a / a² ≥ Δ₀²` at the asymptotically-free scaling rate, i.e. `λ_a ~ a² Δ₀²`, with `λ_a` transported multiplicatively by `Z_O²`. The transport itself is free (it is a square); what is not free is the non-zero seed it multiplies. This distinction pervades the honest reading of the numerics in §5 and §9.

### 2.4 The node board, in one table

| Half | Object | Local completion | Status |
|---|---|---|---|
| G1 | 4D continuum existence | 100% assembly modulo `{L3, G2}` | [T]+[C], one UV node open |
| G1.2 | UV uniform `Z_a` = L3 | ~40% | [M] |
| G2 | Gap: region coverage | ~90% | theorems G2.1–4 |
| G2.5–6 | Gap: sign / proof | **0%** | [M] = H1, the species |
| Channel | `0⁺⁺` identification + no dissolution | 100% | [C] |
| **Clay** | existence ∧ gap | **0%** | **open** |

The percentages are **node-local** — progress toward each node's own goal — and are emphatically **not** a percentage proof of Clay. Clay is 0% because the species is unbuilt.

---

## 3. The exact spectral squeeze, and why H1 is necessary but not sufficient

### 3.1 Setup

Let `ν ≥ 0` be the Källén–Lehmann measure of `tr F²`, with `supp ν ⊂ [m₁², ∞)`, atom weight `f₁² = ν({m₁²})`, and moments `M_k = ∫ s^k dν(s)`. Let `a` be the zero-momentum coarse mode (`a = ⟨v, log U⟩` for the relevant mode; `Var(a)` is the zero-momentum connected propagator). The **susceptibility** is

```
χ := Var(a) = ∫ s⁻¹ dν(s)  (= the zero-momentum connected two-point function) .
```

This is the inverse of the effective curvature only in the Gaussian case; in general it is a residue-weighted *harmonic* mean of the masses, dominated by the lightest residue. It is not the lightest mass.

### 3.2 The squeeze [T]

For any such `ν`:

```
f₁² / χ  ≤  inf supp ν  ≤  M₀ / χ .
```

**Proof of the lower inequality.** `χ m₁² = f₁² + m₁² ∫_{s > m₁²} s⁻¹ dν(s) ≥ f₁²`, since `m₁² s⁻¹ ≥ 1` on the support above `m₁²`. Hence `m₁² ≥ f₁² / χ`. **Proof of the upper inequality.** On `supp ν`, `s ≥ m₁²`, so `s⁻¹ ≤ m₁²⁻¹`, hence `χ = ∫ s⁻¹ dν ≤ M₀ / m₁²`, i.e. `m₁² ≤ M₀ / χ`. ∎

Both inequalities are elementary and exact. They are the entire content of the "squeeze" and they contain no dynamical input.

### 3.3 H1 is necessary, not sufficient

The condition H1 is

```
H1 :   Var_μ(a_{k→0}) = χ ≤ 1/λ ,   uniform in (L, a) .
```

**Necessity.** If the theory has a gap `m₁ ≥ √λ` (with the normalization `M₀` fixed), then `χ = ∫ s⁻¹ dν ≤ M₀ / λ`, so `χ ≤ M₀/λ`. This is a theorem: `gap ⇒ χ ≤ M₀/λ` (up to the normalization constant). The corpus records it as `χ = ∫ s⁻¹ dν ≤ M₀/λ` (`REDUCTION-AUDIT.md`).

**Non-sufficiency.** Take

```
ν = w δ_ε + (1 − w) δ_S ,    w/ε ≤ 1/λ .
```

Then `χ = w/ε + (1−w)/S < ∞`, so H1 holds, while `inf supp ν = ε → 0`. **H1 holds and the gap fails.** The measure has a thin atom at `ε`, invisible to the susceptibility because the susceptibility is a harmonic mean: a tiny residue at tiny mass contributes `w/ε`, which can be kept below `1/λ` by taking `w` small. The susceptibility sees the *residue-weighted* mass, not the lightest mass.

This counterexample is not a pathology; it is the exact obstruction the campaign met everywhere. It is why the audit's verdict is that the chain written as `⟺` is unsound: read against the corpus's own definitions, it is **(iv) false → (iii) conjecture → (ii) definition**, plus one genuine one-directional theorem. The final `⟺` (H1 `⟺` gap) holds in one direction only.

### 3.4 The correct sufficient condition

The repair is the pair:

```
(f₁² ≥ r₀)  ∧  (χ ≤ 1/λ)   ⇒   m₁² ≥ f₁² / χ ≥ r₀ λ > 0 .
```

This is immediate from the squeeze and it is the campaign's central logical content: **a lower bound on the mass needs both a residue floor and a susceptibility ceiling.** The residue floor is the hard half; the susceptibility ceiling is H1. The campaign's honest statement of the missing input is therefore not H1 alone but

```
H*  :   f₁² ≥ r₀ > 0   with isolation I2 ,   (and χ ≤ 1/λ as the matching ceiling),
```

where I2 is `Λ := inf supp(ν − f₁² δ_{m₁²})`, `Λ − m₁² ≥ δ > 0`. The residue–level-gap lemma turns I1 into I2 and back, so `H* = H1` *as a structural identification*; but the identification is the conjecture, and the residue datum is strictly stronger than the bare gap.

### 3.5 Why H1 is blind: Schur concavity

H1 is an **upper-only** condition. Every functional of `ν` built from its moments or from positivity is a Schur-concave functional of the measure, hence gives an **upper** bound on the atom weight, never a lower bound. Concretely, for any nonnegative test function `φ`,

```
f₁² ≤ M₀ ,        f₁² ≤ χ m₁² ,
```

and the first is the trivial "atom ≤ total mass," the second the squeeze. There is no positivity-only route to a *lower* bound on `f₁²`: the free massless field is a legitimate OS theory with `ν` supported on `[0, ∞)` and `f₁² = 0`, and any positivity argument that would force `f₁² > 0` must also force it for the free field, which is false. This is the master obstruction (§6.1), and it is why H1 — a second-moment inequality — is structurally blind to the light state that the squeeze's counterexample exploits.

### 3.6 `G″` versus `1/χ`

The curvature `G″ = V_eff″` and the inverse susceptibility `1/χ` coincide only for a Gaussian marginal. Ward's identity

```
G″ · Var(a) = 1    ⟺    the marginal of a is Gaussian ,
```

is exact in one direction only. Measured at `β = 2.5`, the two differ by 2–5% (`G″ ≈ 0.281` versus `1/Var ≈ 0.272`; `FORGE5-SL-IR-SINGLE.md`, quoted in `docs/new-tools-registry.md`). So even the "curvature face" is not literally the susceptibility; it is a second-moment object that equals the susceptibility only in the Gaussian limit. The unification `G″ = 1/χ = H1` is therefore a *conjecture with a measured 2–5% gap*, not an identity. The exact identity that does hold is the variance decomposition of §4.6.

### 3.7 The verdict on the reduction

The honest logical skeleton (`REDUCTION-AUDIT.md` §"What actually survives") is:

- **[T]** isolated atom at `inf supp ν` of a channel ⇒ exponential clustering in that channel;
- **[T]** `gap ⇒ χ ≤ M₀/λ`;
- **[T]** `G″ ≥ λ ⇒ Var(a) ≤ 1/λ` (Brascamp–Lieb, one way);
- **[C]** the unification `H* = H1` ("one object") — a structural identification, not a theorem;
- **[iv] false** `gap ⟺ atom`, and `H1 ⟺ gap`.

Written as `⟺`, the chain is unsound. Its true shape is a **diamond of implications converging on the necessary condition `χ ≤ 1/λ`**, with the sufficient condition being the strictly stronger `H* = I1 ∧ I2`. Proving H1 alone does not yield the gap (the tiny-residue counterexample); proving H* does, but H* is not H1. This is the first place where the campaign's honesty line is load-bearing, and we state it plainly: **the reduction relocates the residue content; it does not remove it.**

---

## 4. Proven theorems

This section states each theorem with its precise hypothesis and label. Every theorem names the file containing the proof. We separate the unconditional strong-coupling results from the conditional ones and state explicitly, where relevant, that a residue theorem is **not** a gap theorem.

### 4.1 The strong-coupling residue theorem (P3) — [T], regime-limited, I1 only

**Source:** `research/unconditional/PROOF-P3-UNCONDITIONAL.md`; `PROOF-P1P2.md`; `PROOF-P3.md`.

**Theorem (P3).** Fix `a > 0`. There exist `β₁ = min(1/48, c/N) > 0` and, for each `β < β₁`, a constant `r₀(β) = ½ c₁(β) u(β)^p > 0` (`c₁` a positive SU(N) group constant, `u(β)` the fundamental character coefficient, `u ≍ β`), such that for every torus `T⁴_L`,

```
f₁²_L = lim_{t→∞} e^{m₁(L) t} C_L(t) ≥ r₀(β) > 0 ,
```

**uniformly in `L`**.

**Hypotheses and imports.** Wilson action `S = β Σ_p (1 − N⁻¹ Re tr U_p)`, Haar measure, character expansion with `c_r(β) ≥ 0`. The proof imports the Osterwalder–Seiler connected-composite expansion and the Kotecký–Preiss criterion; the single lemma left as "standard assembly" by the earlier drafts — the **graded-strip Kotecký–Preiss lemma** — is written out in full in `PROOF-P3-UNCONDITIONAL.md` §2–§4. Its hypothesis is the per-site condition `u(β) K < 1` with explicit `K ≤ C N` (the plaquette adjacency degree is `≤ 12`; the SU(N) link multiplicity is polynomial in `N`, bounded by `C N` on the fundamental sector). The one constant not reduced to a decimal is the constant in `C_g ≤ C N`; its *existence* is proved, and its value only shifts the prefactor of `β₁`.

**Proof skeleton.** (i) `tr F²` is a finite A1g sum of local plaquette observables, so `C(t)` is a sum over connected polymers; absolute convergence is the KP condition `uK < 1`, uniform in volume and strip length. (ii) A connected polymer linking time `0` to time `t` must meet every intermediate slice, so its plaquette cost is `≥ t/a`; equality forces a single fundamental flux tube, weight `a₀ u^{t/a} = a₀ e^{−m₁ t}`. This tube is the *unique* leading concentration, and its amplitude `a₀(β) > 0` is a product of positive fundamental character coefficients times the positive singlet recoupling `f ⊗ f̄ → 1`. (iii) `f₁² = |⟨Ω|O|G⟩|²` is a modulus-square, so no sign can undo it. (iv) Every non-tube polymer costs a whole factor `u` more; the KP criterion is independent of backbone length, so the remainder is `O(u)` uniformly in `t, V, L`.

**What it is.** A genuine **volume- and L-uniform lower bound on the residue** `f₁²` at strong coupling, fixed `a`. It is the half that Shen–Zhu–Zhu does not supply: SZZ gives a gap (a decay-rate *upper* bound), P3 gives a residue *lower* bound.

**What it is not — the correction that must be preserved.** Earlier drafts claimed the integer cost grading "additionally yields I2 (isolation)." **That claim is wrong** and is corrected in `PROOF-P3-UNCONDITIONAL.md` §5 and `PROOF-P1P2.md` §3. The cost grading grades **polymer amplitudes**, not **transfer-matrix eigenvalue spacing**. It makes the single tube the unique leading term (I1), but the KP/OS remainder is a `t`-uniform `O(u)`, never an `O(e^{−δt})`; so it cannot separate a second atom from one at `ε → 0`, nor place the excited `0⁺⁺*` (at `≈ 1.5 m₁`) below the `2 m₁` continuum. **I2 remains exactly the open step.** Hence P3 is a **residue theorem, not a gap theorem.** The strong-coupling gap itself is classical Osterwalder–Seiler (1978) and is not what P3 proves.

**Scope.** Fixed `a`; `β < 1/48`; no continuum limit `a → 0`; no transport to `β ≈ 2.5`. The analyticity radius `1/48` is the Ursell series radius, not a wall in the spectrum: the residue `f₁²(β) = |a(β)|²` is real-analytic wherever `λ₁` is simple and isolated (Kato; analytic in a disk), so `1/48` is a limitation of the *method*, not of the object. But the method stops there, and no absolute strong-coupling method can cross its analyticity radius.

### 4.2 The no-dissolution kinematic theorem — [C]

**Source:** `research/unconditional/PROOF-NO-DISSOLUTION-KINEMATIC.md`; `PROOF-DISSOLUTION-SPECTRAL-DENSITY.md`.

**Theorem (no dissolution).** Suppose (K-A) the theory has a uniform gap `m_g ≥ μ > 0`, and (K-B) the lightest `0⁺⁺` satisfies `m_{0⁺⁺} < 2 m_lightest` (implied by `0⁺⁺` being the overall lightest). Then the lightest `0⁺⁺` is an isolated atom of its channel, separated from every continuum by a margin equal to its own mass: `m_{0⁺⁺} < θ_{0⁺⁺} = 2 m_g`.

**Proof skeleton.** In the confined phase there is no asymptotic gluon, hence no `(gg)` cut; the lowest continuum in the `0⁺⁺` channel is two glueballs, at `θ ≥ 2 m_g`. If the `0⁺⁺` is the overall lightest, `m_{0⁺⁺} = m_g`, so `m_{0⁺⁺} < 2 m_{0⁺⁺} ≤ θ`. The atom sits strictly below every continuum by its own mass. Dissolution would require a continuum at or below `m_{0⁺⁺}`, i.e. a state lighter than `½ m_{0⁺⁺}` — a level crossing **from below** — or gap closure. Same-channel `0⁺⁺` versus excited `0⁺⁺*` is a von Neumann–Wigner avoided crossing (order cannot flip); a cross-channel `{0⁺⁺, 2⁺⁺}` ordering flip is a relabeling (residues unchanged), and true dissolution needs the factor-2 hierarchy `m_{0⁺⁺} ≥ 2 m_lightest`.

**Sharpening.** In `0⁺⁺`, charge conjugation `C` and parity `P` cannot separate `gg` from `(gg)(gg)` (both are `C`-even, `P`-even), so the *operative* threshold `Λ` for isolation is the excited `0⁺⁺*` at `≈ 1.5 m_g` (physically `1.56 m_{0⁺⁺}`), not the two-glueball `2 m_{0⁺⁺}`. The coarser target "no `0⁺⁺` below `2 m_g`" is false; the correct I2 is a gap to the excited glueball of size `δ ≈ 0.5 m₁`, fixed by power counting, not symmetry (`PROOF-ISOLATION.md` §4; `PROOF-EDGE-E2.md`).

**Conditional.** The theorem is a genuine kinematic theorem but inherits rather than proves its antecedents: (K-A) is the gap, and (K-B) is the ordering input (`IN-11`). It fails only on gap closure (deconfinement, not the zero-temperature smooth crossover) or on a factor-2 hierarchy. Both are dynamics, neither a smooth-crossover effect. The result is therefore [C]: it converts "is the `0⁺⁺` an atom?" into "is it below the two-body threshold?" — and the antecedent is exactly the open dynamics.

### 4.3 Non-negativity of the Haar geometric curvature — [T]

**Source:** `research/unconditional/PROOF-JACOBIAN-SAVES.md`.

**Theorem (W_geo ≥ 0).** Work in the flat chart `X_ℓ = log U_ℓ ∈ su(N)`, so Haar is `dU = J dX` with `Φ = −log J`. For `SU(N)`, `Φ = Σ_{i<j} g(λ_i − λ_j)` with `g(d) = −2 log( sin(d/2)/(d/2) )`; on `|d| < 2π`, `g″ ≥ 1/6`. Hence `Φ` is convex, `e^{−Φ} dX` is log-concave, and by Prékopa every linear marginal `a = ⟨v, log U⟩` has convex free energy:

```
W_geo := E[Φ″] − Var(Φ′) ≥ 0    for every mode, every k including k → 0, uniformly in L.
```

Additionally, pointwise, `Φ″ ⪰ 0` under any measure and any `β`, giving Lemma B: `E[Φ″|0] ≥ (N/12)‖v‖²` (SU(2) `1/6`, SU(3) `1/4`; measured `0.34`, `0.46`). At `β = 0` the geometric term is the whole curvature: SU(2) `V_eff″ = 0.248` versus `E[Φ″] − Var(Φ′) = 0.246`, an exact closure.

**The exact split and what fails.** The identity is

```
G″(0) = E[H_vv + Φ″|a=0] − Var(G_v + Φ′|a=0) = Δ_dyn + W_geo + X ,
Δ_dyn = E[H] − Var(G) ,      X = −2 Cov(G, Φ′) .
```

`W_geo` is the curvature of the pure-Haar marginal; `X` is the cross-fluctuation of the action force with the Jacobian. The theorem bounds `W_geo` from **below by 0**, but the *measured* `W_geo ≈ +0.32` is what carries the sign at `β = 2.5`. Certifying `W > Var − E[H]` needs an upper bound on `Var(G|0)`, which is precisely the induced-mass/H1 object. So this is a **sign theorem, value-free**, not a gap theorem. The `−0.165` "DEAD-lean" value reported elsewhere is a proxy that dropped exactly `+W`; measured exactly on its own mode, the marginal is `+0.51` (§5).

### 4.4 Gauge-invariant Perron–Frobenius — [T] finite volume, plus an edge no-go

**Source:** `research/unconditional/PROOF-GAUSS-INVARIANT-PF.md`.

**Theorem (gauge-invariant PF).** On a finite lattice, `β > 0`, the transfer matrix `T = e^{−aH}` has a strictly positive kernel `T(U,V) > 0` for all `U, V` (a limit of the electric heat kernel, which is positive by `K_t(g,h) = Σ_λ d_λ χ_λ(g⁻¹h) e^{−tC_λ} > 0`, times a positive magnetic weight). Since `T(gU,gV) = T(U,V)`, it descends to a strictly positive kernel `T̄` on `𝒜/G`; for `0 ≤ f ∈ H_phys`, `f ≠ 0`, `(Tf)(U) > 0` everywhere. So `T` is positivity-improving at every `β > 0`, and by Jentzsch/Krein–Rutman `λ₀ = 1` is simple with a strictly positive vacuum `Ω ≡ 1`, and `0 ≤ λ₁ < λ₀`, giving a finite-volume gap `Δ_V = −(1/a) log λ₁ > 0`.

**Why this is honest and why it falls short.** The proof uses the **full positive Haar measure**, no gauge slice, so no Gribov copy or Faddeev–Popov sign enters — the honest content of "a positivity gauge-fixed PF lacked." But (i) `λ₁` is the lightest state over *all* channels, not the `0⁺⁺` specifically; (ii) `inf_V Δ_V > 0` is the mass gap and is not obtained — PF/Dobrushin constants scale like `e^{−cV}`, so finite-volume gaps are *upper* bounds on `Δ_{ℝ⁴}`; (iii) at `a → 0` compactness fails, and the free field (reflection-positive, PF, strictly positive vacuum) has `Δ_V → 0`.

**The exact edge obstruction [N].** The vacuum is `0⁺⁺` and the lightest scalar glueball is *also* `0⁺⁺`: they share the channel. The `0⁺⁺`-projected transfer matrix is still positivity-improving, so its Perron vector is the **vacuum** (top of the sector) and the glueball is the second state, in `Ω^⊥ ∩ (0⁺⁺)` — exactly where the single reflection-positivity cone degenerates. To reach the glueball one would need a second invariant cone / total positivity (a Karlin flag) on `𝒜/G`, but the quotient carries **no order**, so no Gantmacher–Krein/TP cone exists. Gauge invariance reaches the vacuum and the gap, never the glueball. The decisive check is on the physical space itself: free Maxwell is PF and strictly positive with a finite-volume gap, yet its `0⁺⁺` spectrum is a two-photon cut from zero, `f₁² = ν({0}) = 0`. Gauge-invariant PF plus strict positivity plus a finite-volume gap **does not imply** a glueball atom.

### 4.5 The strong-coupling coarse Poincaré lemma (Lemma SC) — [T], regime-limited

**Source:** `research/unconditional/PROOF-COARSE-POINCARE.md` §(c).

**Lemma SC.** At `β < 1/48`, Shen–Zhu–Zhu give a uniform LSI `Var_μ(F) ≤ (1/λ_μ) ∫ |∇F|² dμ`, uniformly in volume. For a block averaging map `P` with `Lip(P) = C_blk = O(1)` (averaging is a contraction, volume-independent), set `F = φ ∘ P`. Then for the pushed-forward measure `ν = P_*μ`,

```
Var_ν(φ) ≤ (C_blk² / λ_μ) ∫ |∇φ|² dν ,
```

uniformly in `(L, a)`. So the coarse Poincaré inequality holds at strong coupling with `λ = C_blk⁻² λ_μ`.

**Scope.** This is a genuine anchor, but it lives exactly where the full LSI already holds. It does not transport to `β = 2.5`, which lies outside every high-temperature regime. The file's three direct routes to `β = 2.5` fail for independent, nameable reasons: the effective action is Savvidy-indefinite (not provably convex; the bare coarse Hessian minimum is `≈ −0.3`, and the one-loop Schur form is undefined because `H_ff` is indefinite); the block gap is a *fast* relaxation rate, not the slow stiffness, and the only variance bound runs the wrong way by `≈ e^{5N_p}`; and the provable regime is `β < 1/48`.

### 4.6 The Feshbach–Bogoliubov two-block lemma and the exact variance decomposition — [T] abstract / identity

**Source:** `research/unconditional/FORGE-BLIND.md` §2; `INVENT-VARIANCE-INEQ.md`; `docs/new-tools-registry.md` A2.

**Lemma (two-block gap; unconditional, abstract).** Let `A = [[A, B], [B*, C]]` on `ℋ = ℋ_- ⊕ ℋ_+`, with `C ≥ γ I`. Define the Schur complement `S := A − B C⁻¹ B*` on `ℋ_-`. If `S ≥ μ I`, then

```
A  ≥  min{ μ / (1 + 2‖B‖²/γ²) ,  γ/2 } · I  > 0 .
```

*Proof.* `⟨u, A u⟩ = ⟨u₁, S u₁⟩ + γ ‖u₂ + γ⁻¹ B* u₁‖²`; put `v₂ = u₂ + γ⁻¹ B* u₁` and use `‖u₂‖² ≤ 2‖v₂‖² + 2‖B‖² γ⁻² ‖u₁‖²`. ∎

The lemma needs **no hypothesis of a positive Hessian**; `A` may be arbitrarily negative. The indefinite Savvidy Hessian is not hidden but *converted* into the demand `S ≥ μ > 0`. For 4D pure SU(N), this becomes the **tachyonic-compensation inequality** (TCI):

```
⟨φ, Q φ⟩ + ⟨φ, K φ⟩ − c F  ≥  μ > 0     for all ‖φ‖ = 1, φ ∈ ℋ_- ,
```

i.e. the gauge-invariant effective potential is convex at its minimum, pointwise and nonperturbatively, and the Savvidy vacuum has no imaginary part. The tool is not vacuous: on the free massless field, `Q = 0`, the bottom of `A|_-` is `0`, `S` is not `≥ μ`, TCI fails, and the certificate correctly returns `0`.

**Lemma NS (eliminating a negative block raises the survivor).** In the quadratic form `[[A,B],[B*,C]]` with `C ≺ 0`, the Schur complement satisfies `Σ = A − B C⁻¹ B* = A + B |C|⁻¹ B* ⪰ A`. This is true and is a genuine gain, but it certifies only the **fast** half of the block BBS scheme; the slow wall is the variance term, which the block elimination never touches (`D-BBS-INTEGRATE-NEG.md`).

**Exact variance decomposition (Savvidy-corrected Brascamp–Lieb).** With `K = Π₀ Hess V Π₀`, `c = Π₀ Hess V v`, `h_n = vᵀ Hess V v`,

```
G″(0) = S(v) − D(v) ,
S = E_μ[ h_n − ⟨K⁻¹ c, c⟩ | 0 ] ,      D = Var_μ(∂_n V | 0) − E_μ[ ⟨K⁻¹ c, c⟩ | 0 ] .
```

Sign law: `D ≤ 0 ⟺ K ⪰ 0` (true BL); `D = 0 ⟺` Gaussian; `D > 0 ⟺` fiber Savvidy. At `β = 0`, `V = Φ` is convex, so `D ≤ 0` and `G″ ≥ S ≥ N/12 > 0` — a theorem. At `β = 2.5`, `D/S ≈ 2%` (measured; the marginal is essentially Gaussian), so `G″ > 0` is numerically robust with ~98% margin. The residual is a uniform bound on the 2% fiber defect, which is H1. The acid test is non-trivial: the free Gaussian has `D = 0`, and the `k → 0` mode has `Var(a) = ∞`, so `S = 0` and `G″ = 0 − 0 = 0`; the tool manufactures no false positivity.

### 4.7 The residue–level-gap lemma and the conditional theorem — [T] / [C]

**Source:** `GAP-G2-edge` (quoted in `G2-ASSEMBLE-100.md` §2); `docs/conditional-theorem.md`.

**Lemma (residue–level-gap).** Let `Λ := inf supp(ν − f₁² δ_{m₁²})`. Then

```
f₁² ≥ (Λ M₀ − M₁)/(Λ − m₁²) = M₀ (Λ − R̄)/(Λ − m₁²) ,      R̄ := M₁/M₀ .
```

The right side is positive **iff** `Λ > R̄`, and `Λ > R̄` suffices (if `f₁² = 0` then `Λ = m₁²` and `R̄ ≥ Λ`). The converse fails (a positive atom with a heavy tail can have `R̄ > Λ`). So this is a *self-certifying sufficient certificate*, not an equivalence. By this lemma and Theorem A of `mass-spectrum`, `atom(I1) ∧ I2 ⟺` sharp two-exponential decay `C(t) = f₁² e^{−m₁t}(1 + η(t))`, `|η(t)| ≤ B e^{−δt}`.

**The conditional theorem.** Under the twelve named inputs `IN-0 … IN-11` (`docs/conditional-theorem.md`), there exists a gauge-invariant measure `μ` on `ℝ⁴` satisfying OS0–OS4, non-Gaussian with `g_ren/Λ > 0`, with `Δ/Λ ≥ min(δ_W, δ_S) > 0`, unique vacuum, area law, and group-uniform constants; and the gap is equivalent to `inf spec(H_a|_{Ω^⊥}) ≥ δ Λ` uniformly in `a`. The synthesis chain is **9/9 closed** (each gluing lemma is either proved or reduced to a named input), and the core of S3, S4, S5, S10, S11 is machine-verified in Lean with frozen axioms (`Handoff.lean`, `RGHandoff.lean`, `HandoffMargin.lean`, `GlueA.spectral_decay`, `GL2a.lean`, `GL2bPrime.lean`, `GL7Dual.lean`, `GroupUniform.lean`, `GL6Fubini.lean`, `Conditional.lean`; the checks are listed in `conditional-theorem.md` §6).

**Honest status of the inputs.** Twelve inputs: two proved or machine-verified (`IN-1` bounded class, vacuous, and `IN-8` assembly), seven open (`IN-2`, `IN-3`, the 4D part of `IN-4`, `IN-5`, `IN-6`, the quantitative constants of `IN-7`, `IN-9`), one explicit magnitude bridge (`IN-10`, `hcert: cert ≤ Δ/Λ`), and one physical input (`IN-11`, the lightest bulk state in `0⁺⁺/2⁺⁺`). `IN-3` — the nonperturbative beta-function bound in the intermediate-coupling region `λ ∈ [1, 7.3]` — is **the single physical wall**. Balaban-type RG stops at `λ ~ 1`; strong-coupling expansions require `λ` large; there is no rigorous method in between. The theorem is [C], and the conditional reads as: *if* these twelve inputs hold, Clay follows — not that they do.

**A crucial correction (red team).** An earlier campaign summary overstated "the certificate is the gap." The correct statement (`campaign-final.md` §6) is that the certificate is a **sufficient** route (`cert + IN-10 + IN-11 ⇒ gap`), but `cert > 0` is strictly weaker than the gap, and the necessity `gap ⇒ cert` is unproved. So "stuck at the certificate" does **not** imply "stuck at the gap"; routes to the gap outside the certificate family are not ruled out. This paper preserves that correction.

### 4.8 Summary table of proved theorems

| Theorem | Label | Hypothesis | What it gives | What it does *not* give |
|---|---|---|---|---|
| P3 residue (I1) | [T] (+OS/KP imports) | fixed `a`, `β < 1/48`, `T⁴_L` | `f₁² ≥ r₀(β) > 0`, uniform in `V, L` | I2 isolation; `a → 0`; `β ≈ 2.5` |
| No dissolution | [C] | uniform gap + `0⁺⁺` lightest | `m_{0⁺⁺} < θ = 2m_g`; atom pinned | the gap + ordering themselves |
| `W_geo ≥ 0` | [T] | Haar, flat chart, `g″ ≥ 1/6` | sign of the geometric term, all k, all L | dominance over Savvidy at 2.5 |
| `E[Φ″] ≥ N/12` | [T] | pointwise, any measure | floor on the Jacobian curvature | strict gap |
| Gauge-invariant PF | [T] finite V | `β > 0`, positive Haar pushforward | vacuum, finite-volume gap | first excited `0⁺⁺` (edge) |
| Lemma SC | [T] | `β < 1/48`, SZZ LSI, 1-Lipschitz blocks | coarse Poincaré, `(L,a)`-uniform | transport to 2.5 |
| Feshbach–Bogoliubov | [T] abstract | two-block, `C ≥ γI`, `S ≥ μI` | explicit gap, no positive Hessian needed | TCI (the 4D input) |
| `G″ = S − D` | [T] identity | exact | `β = 0` positivity `G″ ≥ N/12` | uniform bound on `D` at 2.5 (= H1) |
| Residue–level-gap | [T] | moments available | self-certifying sufficiency `Λ > R̄` | equivalence; `Λ` lower bound |
| Conditional theorem | [C] | `IN-0 … IN-11` | Clay conclusion | the inputs (esp. `IN-3`) |

---

## 5. The six tools and their free-field acid tests

The project charter imposes a **free-field acid test** on every purported new tool: it must fail on the free massless field, where `f₁² = 0` and `G″ = 0`. A tool that "proves" a gap for the free field is vacuous and must be discarded. Six constructions pass.

### 5.1 Multiplicative residue transport

**Formula.** `f₁²(μ) = Z_O(μ, μ₀)² · f₁²(μ₀)`, with `Z_O = (β/2g)(μ₀)/(β/2g)(μ) > 0` and anomalous dimension `γ_O = −d log(β/2g)/d log μ = 2 b₀ g² + …`.

**Basis.** `O_S = ¼ tr F²` renormalizes multiplicatively: the identity mixes diagonally, `F F̃` is `P`-odd, and BRST-exact derivatives vanish on physical states.

**Acid test.** Passes: `0 → 0` (multiplicative); the free `:φ²:` has `f₁² = 0` by `Z₂` parity (`a ≡ 0`); the constant `Z_O² > 0` because it is a square, and is `O(1)` in `V`.

**Residual.** Conditional on `(CP)`: a single `T = 0` confining phase, a `0⁺⁺` edge that is Kato-preserved, and uniform isolation. The theorem `f₁²(2.5) ≥ Z_O² r₀(β_eff)` is **[C]** on `(CP)`. It transports the *sign and positivity* of the residue; it transports neither the size nor the isolation.

### 5.2 Exact variance decomposition (Savvidy-corrected Brascamp–Lieb)

**Formula and acid test.** Stated in §4.6. The decomposition `G″ = S − D` passes the acid test in the strong sense: the free Gaussian has `D = 0`, the `k → 0` mode has `Var(a) = ∞ ⇒ S = 0`, and `G″ = 0 − 0 = 0`; no false positivity is manufactured.

**Residual.** A uniform bound on the fiber defect `D` at `β = 2.5`; measured `D/S ≈ 2%`. This is H1 in curvature language.

### 5.3 The Fredholm determinant identity

**Formula.** `f₁² = −λ₁ Res_{z = z₁} G(z) = −(1/z₁) N_O(z₁)/D_O′(z₁)`, with `D_O(z) = det(1 − z T_O)`, `log D_O = −Σ_k z^k Tr(T_O^k)/k`, `z₁ = e^{a m₁}` the smallest-modulus zero, and "periodic orbits" `Tr T_O^k = Z_{S¹_k × T³}` (closed-time flux tubes).

**Acid test (strict).** The free field has a pure cut and no pole, so its residue is `0`; YM has an isolated zero below its own two-glueball cut (`m₁² < 4 m₁²`), so `f₁² > 0`.

**What it buys.** The pole and the cut have identical lower moments; only the zero set differs. So the determinant is **linear-blind** but **sees the atom**, and it supplies `m₁²` (the zero location), which no moment does.

**Residual.** Evaluating the zero at `a → 0` is the marginal BPHZ wall; the continuum determinant's renormalization *is* the marginal gauge BPHZ, and the zero coalesces into the Källén–Lehmann pole (`CONTINUUM-FREDHOLM.md`).

### 5.4 Gauge-invariant IR mass

**Formula.** `m_GI := inf supp ρ_s = lim_{|x|→∞} −∂_{|x|} log|G_s(x)| = E₁^{(s)} − E₀`, on the physical Hilbert space `H_phys = L²(𝒜/𝒢, μ̄)`, with `s = (β/2g) tr F_{μν}F^{μν} ∈ H_phys` and `ρ_s ≥ 0` by reflection positivity.

**What it removes.** The claim that the gauge potential `A ∉ H_phys` and that `ρ_A` has indefinite sign are both artifacts of using `A`; with `F²` the spectral density is non-negative. **Gauge invariance was never the wall.**

**Acid test (definitional).** Free Maxwell has `m_GI = 0` (a threshold two-photon cut from zero); YM has `m_GI = m_{0⁺⁺}`.

**Residual.** The edge atom `f₁² > 0`; the definition removes a pseudo-obstruction but does not manufacture the atom.

### 5.5 Gauge W-entropy

**Formula.** `F_W[μ_t] = ∫_𝒜 (S_W(A) + |∇_A f_t(A)|²_G) e^{−f_t(A)} dvol_G`, along the conjugate flow `∂_t A = −∇_A S_W`, `∂_t f = −Δ_A f + |∇_A f|² − S_W`. Gap preservation is monotone iff `Γ₂ = Ric_𝒜 + Hess_𝒜 S_W ≥ ρ`.

**Acid test (mirror).** Free field: `Γ₂ ≥ 0`, monotone, but it certifies only the free lattice gap (non-empty). YM at `β = 2.5`: `Γ₂ ≈ −2.0 < 0` — correct on the gapless side, **blocking** on the gapped side.

**Verdict.** Perelman's `W` requires "geometry = flow" to be the same object; in flat-space YM they split. What must be crossed is not a new functional but the missing non-perturbative flow itself.

### 5.6 Feshbach / Krein

The Feshbach–Bogoliubov certificate is §4.6: an unconditional abstract two-block lemma whose only debt is TCI. The **Krein-space** analysis is the precise language for the indefinite Savvidy Hessian `−D² + 2 ad F_B`: Pontryagin spaces, Langer definitizable operators, spectral functions. It rigorously proves the UV tower cone-empty statement (shells `λ_j > k_⋆` have `κ = 0`, `Tr H_- ≲ ∫ |F|³`) — a genuine **UV half-theorem**. But every Krein invariant is a **count/sign** (negative index `r_-`, signature) or a **UV scale** (the `J`-positive floor `= Θ(g|B)`, Landau levels), whereas the gap is the **min** at the bottom of `k → 0`, an **IR** scale. The two are mutually exclusive: choosing `J = sign H` makes the floor physical but gives `κ = dim H_+ = ∞` and is circular (needs a positive spectrum); choosing an a priori BRST/Gupta–Bleuler metric gives `κ = O(1)` but is literally identical on the free field (acid-vacuous). A signed negative index is an integer, hence rigid (B3), hence equal to its free value `0`. The definitizing split loses uniform margin exactly at `k = 0` (the Savvidy band bottom, where the coarse mode and the most negative mode coincide) — the Krein critical point is precisely the point at which the gap must be read. The residual is `inf spec(H|_{physical side}) ≥ κ_IR > 0`, `(L,a)`-uniform, which is H1.

### 5.7 Scorecard

| Tool | Formula / object | Acid test | Touches the atom? | Exact residual |
|---|---|---|---|---|
| Multiplicative transport | `f₁²(μ) = Z_O² f₁²(μ₀)` | pass (`0→0`, square) | sign/size transport | `(CP)` |
| Variance decomposition | `G″ = S − D` | pass (free `0−0`) | curvature sign | uniform `D` at 2.5 = H1 |
| Fredholm determinant | `f₁² = −λ₁ Res G` | strict pass | **yes** (zero set) | `a → 0` zero |
| Gauge-invariant IR mass | `m_GI = inf supp ρ_s` | pass (Maxwell 0) | removes artifact | edge atom |
| Gauge W-entropy | `F_W`, `Γ₂ ≥ ρ` | mirror pass | no (gap side blocks) | missing flow |
| Feshbach/Krein | two-block lemma, `Π_κ` | pass (free returns 0) | reduction to TCI | TCI = H1 |

Three of the six — transport, variance decomposition, determinant — genuinely "touch" the atom; all three leave the wall on the far side of "continuum plus transport to 2.5," which is the analyticity/BPHZ wall.

---

## 6. No-gos: what each class cannot do

A no-go is a theorem that a *class* of approaches cannot succeed, not a statement that no approach can. The campaign proved several, and they are as valuable as the theorems because they close branches.

### 6.1 Schur-convex master obstruction — [N]

**Statement.** Any `g`-continuous, protected functional of the measure that vanishes at `g = 0` (the free field) cannot produce a positive atom at `g > 0` by positivity or analyticity alone. The free field satisfies all OS axioms and has `f₁² = 0`; any positivity-only argument forcing `f₁² > 0` would force it there too.

**Consequence.** Positivity, reflection positivity, spectral positivity, and analyticity give **upper bounds and the vacuum**, never the edge. This is the single most important no-go; it is why the campaign always says "positivity is blind to the edge."

### 6.2 Goldstone uniqueness — [N]

The only positivity-plus-analyticity theorem that forces an atom is Goldstone's theorem, and it lives at **zero mass** (broken continuous symmetry). The atom `f₁²` is the "positive-mass cousin" of a Goldstone mode, and the positive-mass cousin does not exist at the level of a general theorem. The one protected positive-mass atom is BPS, but it is in the **wrong channel** (central charge/flux; `Z = 0` on `tr F²`).

### 6.3 Periodic no-go (no convexification) — [N]

**Statement.** The log-density of `μ` is periodic on the compact group; a periodic concave function is constant, so there is no log-concave substitute equal to `μ`. The Savvidy negativity is the cost of periodic winding at `|x| > π/2`. Thermal smoothing gives only semi-concavity, and the Prékopa lift is blocked by the indefinite `G″/H_ff`.

**Consequence.** Convexification of the effective action is impossible; the wall is real, not a gauge artifact. (W4 measured 7.1% negative directions after the gauge quotient.)

### 6.4 Nested-cone no-go — [N]

**Statement.** Every sign-changing / Karlin cone `K₀ ⊂ K₁ ⊂ …` contains the positive cone, so on each layer it swallows the strictly positive free field and is vacuous. A cone that passes must cut *below* positivity (where the symbol is indefinite), and there it presupposes the very atom it is meant to prove (circular).

**Consequence.** The "second cone" route to the edge is closed. Plaquette products do not tensorize to a sign-regular cone.

### 6.5 Rigidity no-go (no index) — [N]

**Statement.** Any index that is a discrete integer and locally constant under symmetry-preserving deformations is constant along the connected `β` family (which reaches the free Gaussian at `β = ∞`), hence equals its free value, hence is `0` for YM. The gap must be a **`β`-non-analytic jump, not any discrete invariant.**

**Consequence.** Noncommutative geometry, localization, BPS, and Witten-spinor indices all return `0` at the free point; none can carry the gap. (The only acid-passing Gribov–Faddeev–Popov index counts gauge copies, not physics.)

### 6.6 Feshbach / Schur angle no-go — [N]

**Statement.** Block-decompose `H = [[A,B],[Bᵀ,C]]` with `C ≺ 0`. Then `Σ = A − B C⁻¹ Bᵀ = A + B|C|⁻¹Bᵀ ⪰ A`: eliminating a negative block **only adds** to the survivor (Haynsworth inertia; numerically verified `λ_min(Σ − A) = −6·10⁻¹⁶`). So "a small angle preserves `Σ > 0`" is a sign error: `Σ > 0 ⟺ A = PᵀHP > 0`, and `A` *is* the gap. The only way an angle enters is through the overlap `ε = sin²θ`, requiring `ε < ½` (`θ < 45°`); on a frozen constant background `θ = 0` but the statement is vacuous (`F ≡ 0`, `g`-dependent), and on the physical `k = 0` mode `ε = O(1)` (measured `PᵀHP` minimum `−0.3`) with no momentum separation. The whole route is `⇔ ε < ½ ⇔ H1`.

### 6.7 Krein / Pontryagin no-go — [N]

**Statement.** Krein counting/sign invariants live at the UV scale or count signs; the gap is the IR **min**. The two are mutually exclusive, as detailed in §5.6. The residual is H1.

### 6.8 Wrong-half, wrong-direction, wrong-regime — [N] (umbrella)

Three structural failure modes recur across families and are worth naming as no-gos of method:

- **Wrong half (UV ⊥ gap).** L3, H2-global, gradient flow, Dirichlet comparison, concentration (Herbst), heat kernel, validated numerics, step scaling, effective string, `S³` — these control the UV, a fast rate, or a finite box, but never the IR spectral edge.
- **Wrong direction (upper bounds).** Bogoliubov, Hylleraas, finite-size scaling, Temple/Weinstein, GEVP — variational methods give **upper** bounds on eigenvalues; a lower bound needs a threshold oracle, and the oracle is circular.
- **Wrong regime (strong coupling only).** Rigorous large-N (Chatterjee), SZZ, random surfaces/CPS, center vortices — these stop at `β < β₀` or `λ < 1/24` and cannot reach `β ≈ 2.5`.

### 6.9 Summary of no-gos

| No-go | Class it kills | Precise reason |
|---|---|---|
| Schur-convex | positivity/analyticity-only | free field has `f₁² = 0`; only upper bounds |
| Goldstone | positivity-forced atoms | the only such theorem is at zero mass; BPS wrong channel |
| Periodic | convexification | periodic concave = constant; Savvidy is real |
| Nested cone | second-cone construction | every cone contains the positive cone → swallows free field |
| Rigidity | index/topological | discrete locally-constant index = free value = 0 |
| Feshbach angle | block-elimination | `Σ ⪰ A`; `Σ > 0 ⇔ A > 0` = gap |
| Krein | UV spectral invariants | count/sign ≠ IR min; split fails at `k = 0` |
| Wrong-half/direction/regime | method families | control the wrong scale, sign, or coupling region |

---

## 7. The one-species reduction and the immune map

### 7.1 The species

Every honest route terminates at the same inequality. Written in its most general form (the `(◆)` of `PROOF-JACOBIAN-SAVES.md`), the missing input is a **single-sided upper bound on a conditioned variance**:

```
(M)  Var_μ( ⟨v, ∇S⟩ | 0 )  ≤  E_μ[ H_vv | 0 ] + W − δ ,      uniform in (L, U′) ,
```

equivalently `Var(a_{k→0}) ≤ 1/λ` (H1), equivalently `m_eff²(k→0) ≥ λ`, equivalently `G″(0) ≥ λ`. The **variance-free face** (the only one that could pierce H1) is a deterministic **negative-association / anti-FKG** inequality

```
Cov_μ( G, Φ′ | 0 )  ≤  0
```

for the continuous, frustrated (cos-sum Wilson action), non-attractive SU(N) Haar × Wilson measure. The campaign proves this object is **absent from every known dependence class** (`C5-NEG-ASSOC.md`, `C6-DEPENDENCE-CLASS.md`): negative association (NA), negative regression/NOD, strong Rayleigh, log-submodular, determinantal (DPP), total positivity (TP₂), Chebyshev/Fréchet, FKG/Griffiths, Brascamp–Lieb, and Simon–Lieb all fail, each for a specific structural reason (no order on the quotient, no real-stable polynomial, no DPP kernel, no common scalar, FKG wrong-signed, BL needs `Hess ⪰ 0` = the Savvidy wall, Simon–Lieb is sign-blind variance = H1, Stein–Hellfer–Sjöstrand constant is `1/λ` = H1). The measured `Cov ≈ −0.11` at `β = 2.5` is the chart-local log-cut shadow of *proven* `Φ`-convexity, not a law: the per-mode `X ∈ [−0.09, +0.16]` changes sign, and `β = 5.7` flips it.

### 7.2 Two scales, one shape

The species appears twice:

- **UV face (L3).** `sup_a Z_a⁻¹ ‖H_a‖_D < ∞`: a single-sided force-variance bound on the `Λ_UV` tower. This is G1's open node.
- **IR face (H1).** `Var_μ(a_{k→0}) ≤ 1/λ`, `(L, a)`-uniform. This is G2's open cell.

The shapes coincide because `Var ≥ 0` can only supply a lower edge, while in `G″ = E[H] − Var + W` the variance subtracts. A tool that upper-bounds a conditioned variance at one scale upper-bounds it at the other; the campaign states that building it would close L3 and H1 **simultaneously**, hence Clay. This is the "one species, two scales" reduction. The unification is the conjecture; the dictionaries (Ward, residue–level-gap, the identity `G″ = E[H] − Var + W`) are theorems.

### 7.3 Why the species is not a repackaging

The project charter explicitly forbids repackaging an old inequality. The species is not: it is neither a positivity statement (those give upper bounds on atoms, §6.1), nor an FKG-type positive association (wrong sign), nor an LSI/BL (needs a convex potential, false at Savvidy), nor a moment/Carathéodory bound (Schur-concave, gives atom upper bounds), nor a sign-blind variance identity (Simon–Lieb). It is a **deterministic negative-association inequality for a continuous, frustrated, non-attractive measure** — by the campaign's own audit, a dependence theorem that is **unknown to exist**. This is the "sixth species" the charter calls for: the first five paradigms (positivity, convexity, spectral theory, Perron–Frobenius, RG) are exhausted, and the sixth is still unbuilt.

### 7.4 The immune map

The campaign's archive holds 73 batch reports across the families `ABC*`, `C-*`, `D-BBS-*`, `LOC-*`, `MS-*`, `LOWER-*`, `FORGE-*` (and `FORGE5-*`), `K0-*`, `BS-*`, `TRANSMUTE-*`; several contain multi-route audits (the `LOWER` report alone checks eight mechanisms), so the distinct-route count is ≈92. The earlier master immune map counted ~60 (`YANG-MILLS-MASTER.md` §5); the batch waves added the rest. **Every one terminates at the same named species.** The table below gives each root family and the precise place it dies.

| # | Root family (representative files) | What it tries | Where it dies |
|---|---|---|---|
| 1 | Free-field-blind positivity / analyticity (`INEQ-T1..T5`, `INEQ-R1..R10`, Pólya, Nevanlinna, Bogoliubov) | force an atom from positivity | `g = 0` free field has `f₁² = 0`; Schur-convex master obstruction; upper bounds only |
| 2 | Cross-channel (`BPS`/SUSY, 1-form anomaly, noninvertible symmetry, dual Meissner) | protection from a symmetry | wrong channel: line/flux/charged sector, or spin-1; local `tr F²` is blind; BPS has `Z = 0` |
| 3 | Non-abelian sign (random flow, random cluster/FK, FKG, spectral independence, Picard–Lefschetz thimbles) | positive/signed cluster calculus | no lattice order on the quotient; SU(N) Weingarten face weights signed; FKG wrong-signed |
| 4 | Wrong half, UV (`L3`, `H2`-global, gradient flow, Dirichlet comparison, concentration, heat kernel, validated numerics, step scaling, effective string, `S³`) | control the ultraviolet or a fast rate | orthogonal to the IR spectral edge; finite box = upper bound; wrong kernel |
| 5 | Rename = gap (`H1`, NEO/ESI, `I2`, Hamiltonian binding energy, one-loop `L3`) | identify with a solvable object | the object is the gap restated, in another vocabulary |
| 6 | Upper-bound direction (Bogoliubov, Hylleraas, FSS, Temple/Weinstein, GEVP) | a variational bound | variational methods give upper bounds; threshold oracle circular |
| 7 | Strong-coupling-only (rigorous large-N / Chatterjee, SZZ, random surface/CPS `λ < 1/24`, center vortex) | extend a proven regime | cannot reach `β ≈ 2.5`; method radius is the obstruction |
| 8 | Index = 0 (NCG, localization, BPS, Witten spinor) | a topological invariant | rigidity no-go: discrete locally-constant index = free value 0 |
| 9 | Exact reduction / duality (`ABC` A-leg: large-N master field, Makeenko–Migdal, 2D Migdal, integrable, `N = ∞ → 3`, Eguchi–Kawai, holography, stochastic quantization) | solve an equivalent model | wrong channel, wrong direction (`S_N ≥ S_∞` reversed), or wrong regime; all eight → H1 or misdirection |
| 10 | Constructive coercivity (`ABC` B-leg) | partition all of config space | region-complete but sign-null: open cell has all four handles null; `G″` = conditioning defect of equipartition; CS `|X| ≤ 0.40 > W_geo ≈ 0.32`, short by `O(0.1) = λ` |
| 11 | New inequality / capacity (`ABC` C-leg, `C-BBS-*`) | a convex floor | WSTC/STC close the cones; open cell needs negative association, proven absent; capacity floor is the Poincaré constant of the 1D marginal `= 1/χ = G″ = H1`, necessary not sufficient |
| 12 | Fixed-cell localization (`LOC-*`) | a cell-local conductance floor | `w ~ ξ` circularity; free massless satisfies every premise; single-cell control does not tensorize to `L → ∞` |
| 13 | Transport / monotone (`MS-*`) | carry the strong-coupling gap across Savvidy | every transport is one-sided (upper bounds) or local (a `β`-disk); the missing inequality is `(SGL) m⁴ ≥ κ⟨¼ tr F²⟩`; residue positivity ≠ gap positivity |
| 14 | Lower-bound machines (`LOWER-AMPLITUDE-MECHANISMS`: GKS, FKG, Simon–Lieb, reverse-BL, Carathéodory, RP-cone, Ginibre, Gaussian correlation) | lower-bound the amplitude | "lower-bound machine, upper-bound target": a lower bound on correlation = lighter mode = **anti-gap**; the only correct lower bound (residue weight) is Schur-concave-forbidden |
| 15 | Feshbach / Krein / squeeze / Dirac-eta / spectral flow (`FORGE-*`) | eliminate the negative block | two-block lemma reduces to TCI = H1; Krein count ≠ IR min; spectral-flow/eta index rigid; squeeze negative block only helps the survivor |
| 16 | `k = 0` effective curvature (`K0-*`) | strict positivity of the coarse curvature | bare Savvidy negative is a **saddle**, not the vacuum; finite-`V` strict positivity is trivial; the constant mode carries toron flat directions; uniform strictness = gap |
| 17 | Bethe–Salpeter two-gluon binding (`BS-*`) | solve the bound-state equation | kernel sign-indefinite (Savvidy); Temple needs the threshold oracle (circular); RP sign-blind; area law is the wrong sector |
| 18 | String-tension → glueball transmutation (`TRANSMUTE-*`) | `σ > 0 ⇒ m_{0⁺⁺} ≥ c√σ` | sector (charged 1-form asymptotic vs neutral 0-form short-distance) and scale mismatch; the available inequality points the wrong way |

**Meta-conclusion.** This is not "we failed to find a route." It is roughly ninety independent frontal attacks, each of which proved the same irreducibility from a different direction. That is near theorem-grade *exhaustiveness evidence* — not a theorem, and not a proof, but a structural fact about where the difficulty lives.

### 7.5 The newest six-route waves

The most recent waves mount several routes each (TRANSMUTE 3, BS 6, K0 6, FORGE 6, LOWER 6, MS 6), all landing on the species:

- **`TRANSMUTE`** (string tension → glueball): the 1-form-to-0-form bridge fails; sector/scale mismatch plus wrong-direction inequality; no non-circular `σ ⇒ m` map.
- **`BS`** (two-gluon Bethe–Salpeter binding): the composite gives no new lever; kernel sign-indefinite; every route ends at `m_g > 0` (H1) or an excited-glueball isolation (I2).
- **`K0`** (`k = 0` positive curvature `A > 0`): object-identification correction — the bare Savvidy negativity is a saddle, not the vacuum — and all routes return to H1.
- **`FORGE`** (forge the missing species): all six wall at `A(k=0) > 0 = H1`; the two-block lemma is the best reduction but its debt is TCI.
- **`LOWER`** (lower-bound the edge / one-sidedness): all six wall at `(MS)`; the lower-bound direction is anti-gap.
- **`MS`** (transport across the crossover): all six reduce to `(SGL)` or an equivalent H1 statement.

---

## 8. What the missing tool must satisfy

The campaign's sharpest deliverable is a specification of the object it could not build. The species must be **one** tool satisfying **five** conditions, stated verbatim from `docs/PROOF-SKELETON.md` §3 and `docs/new-tools-registry.md`:

1. **Nonperturbative.** It must be valid outside every convergent expansion, at `β ≈ 2.5`, where the relevant `g ≈ 1.5` and no small parameter exists. (A perturbative resummation is excluded by the Savvidy negative mode and the IR renormalon/Borel ambiguity of the `a⁻⁴` coefficient.)
2. **Sign-controlling lower bound.** It must produce a *lower* bound `m_eff²(k→0) ≥ λ > 0` (equivalently `G″ ≥ λ`, `f₁² ≥ r₀`), not another upper bound. This is the rarest direction: almost every physical tool gives upper bounds or sign constraints.
3. **Non-circular.** It must not assume the gap, the threshold, or the isolation it is meant to prove. (Temple/Weinstein, the BS threshold, the second-level `E₂`, and the cutoff-independent Poincaré floor all fail here by presupposing the answer.)
4. **Non-SUSY and fails the free field.** It must not rely on supersymmetry (the physical theory is pure YM) and it must **fail** on the free massless field, where `f₁² = 0`, `G″ = 0`. This is the acid test: a tool that "proves" a gap for the free field is vacuous.
5. **Uniform lower bound.** It must be uniform in `(L, a)` — volume- and cutoff-uniform — and survive the continuum limit. A fixed-`a` or fixed-`L` statement does not give Clay.

### 8.1 The equivalent forms of the object

These five conditions define the same object in several vocabularies (`docs/YANG-MILLS-MASTER.md` §6):

1. **H1:** coarse Poincaré `Var_ν(φ) ≤ (1/λ) ∫ |∇φ|²`, uniform `λ > 0`.
2. **(CP):** a single `T = 0` confining phase (no bulk transition) + a `0⁺⁺` edge that is Kato-preserved + uniform isolation.
3. **Edge atom `f₁² > 0`:** `ν({m₁²}) > 0` (usc, non-diagonal, needs `E₂`).
4. **`G″(0) > 0` at `k → 0`:** the coarse curvature, `= S − D`, with `D/S ≈ 2%`.
5. **Determinant zero in the continuum:** `det(1 − z T_O)` at `a → 0`.
6. **Savvidy sign** (non-abelian): the indefinite Hessian's sign resolved.

These are not six problems; they are one — the volume-uniform lower bound on `m_{0⁺⁺}` along the path, which is the mass gap itself. The campaign's final convergence (`docs/PROOF-SKELETON.md`) further reduces them to **one cell**: the "typical field × high momentum" cell of the `(field × momentum)` plane, requiring `Var_μ(⟨v,∇S⟩|0) ≤ W_geo`. Every other cell is closed: low field by `Φ″ ⪰ N/12` + BL, the high-field cone by WSTC (measure starvation), the AF-UV cone by `g → 0`. **That cell is a variance-upper-bound species, the same species as G1's L3, at a different scale.**

### 8.2 The machine that does not exist

To close Clay one must build **a nonperturbative machine that transports the strong-coupling gap and residue continuously across the Savvidy crossing and takes `a → 0`**. The campaign states the verdict plainly (`docs/yang-mills-frontier.md` §7, `YANG-MILLS-MASTER.md` §9):

> The last step is to give a volume-uniform lower bound on `m_{0⁺⁺}` across `[1/48, 2.5]`. Every tried tool stops at the strong-coupling end (SZZ, large-N) or on the wrong half (UV). What is needed is a machine that carries the strong-coupling gap across the Savvidy/sign crossover continuously to the physical point — **this machine does not exist.** It is equivalent to nonperturbative 4D dynamics, the classification of nonperturbative pure-gauge 4D CFTs (`TOOL-not-iii.md`).

The strongest contemporary rigorous frontier (CPS/CNS 2025) reaches only `λ < 1/24` (the older "`75×`" figure was a unit confusion; see `PROOF-EMERGENT-4D.md`), nowhere near `β ≈ 2.5`. The machine is the sixth species, and it is unbuilt.

### 8.3 What building it would close

A single one-sided variance upper bound, nonperturbative, coupling-free, `(L,a)`-uniform, failing at `f₁² = 0`, would close **both** G1's L3 (UV) and G2's H1 (IR) simultaneously. That is the precise sense in which "one species, two scales" is the whole gap between the campaign and Clay.

---

## 9. Honest verdict

### 9.1 The headline

**Clay is 0% — open.** This is not a proof of the Yang–Mills mass gap, and no step of this paper should be read as one. The Clay problem is a Millennium Prize problem and remains unsolved.

What is true:

- **G2's region coverage is ~90%**, closed by genuine theorems (G2.1–G2.4); **its sign is 0%**, the open cell being exactly H1.
- **G1 is 100% assembled modulo `{L3, G2}`**, with the one genuinely open UV node L3 at ~40%, and a coarse aggregate of ~70%.
- **The channel half is 100%** (conditionally): the no-dissolution theorem pins the lightest `0⁺⁺` below its continuum.
- **The reduction is exact and honest:** the mass is squeezed as `f₁²/χ ≤ inf supp ν ≤ M₀/χ`; the correct sufficient condition is `(f₁² ≥ r₀) ∧ (χ ≤ 1/λ)`; H1 alone is necessary but not sufficient.

### 9.2 What was actually banked

Positive ledger (all traced to files):

1. **P3 strong-coupling residue theorem** (`PROOF-P3-UNCONDITIONAL.md`): `f₁² ≥ r₀(β) > 0` for `β < 1/48`, volume- and L-uniform. Label [T] modulo standard OS/KP imports. **I1 only; not I2.**
2. **No-dissolution kinematic theorem** (`PROOF-NO-DISSOLUTION-KINEMATIC.md`): the lightest `0⁺⁺` is pinned below every continuum by its own mass. Label [C] on the gap + ordering.
3. **`W_geo ≥ 0`** (`PROOF-JACOBIAN-SAVES.md`): the Haar geometric curvature is non-negative on every mode including `k → 0`, every L. Label [T].
4. **Gauge-invariant Perron–Frobenius** (`PROOF-GAUSS-INVARIANT-PF.md`): finite-volume vacuum and gap, no Gribov/FP sign. Label [T] finite volume; edge obstruction [N].
5. **Lemma SC** (`PROOF-COARSE-POINCARE.md`): strong-coupling coarse Poincaré, `(L,a)`-uniform. Label [T], regime-limited.
6. **Feshbach–Bogoliubov two-block lemma** (`FORGE-BLIND.md`): explicit gap with no positive-Hessian hypothesis. Label [T] abstract.
7. **Exact variance decomposition `G″ = S − D`** (`INVENT-VARIANCE-INEQ.md`): `β = 0` positivity; `D/S ≈ 2%` at `β = 2.5`. Label [T] identity + [M] at 2.5.
8. **Residue–level-gap lemma** (`GAP-G2-edge`): a self-certifying sufficient certificate. Label [T].
9. **Conditional theorem** (`docs/conditional-theorem.md`): `IN-0 … IN-11 ⇒ Clay`, 9/9 gluing lemmas closed, partial Lean verification. Label [C].

Negative ledger (the no-gos): Schur-convex, Goldstone, periodic, nested-cone, rigidity, Feshbach-angle, Krein, and the wrong-half/direction/regime umbrella. Plus three campaign-specific no-gos: **periodic** (no log-convexification), **nested** (cone swallows the free field), **rigidity** (index = 0; the gap is a `β`-nonanalytic jump).

Tools: six constructions pass the free-field acid test, three of them genuinely touch the atom, all three leave the wall on the far side of the continuum limit.

### 9.3 The numerics, read honestly

- **`f₁² > 0` is seen where it is resolvable.** At `β = 2.2` on SU(2), `L = 6`, `N = 4000`, the connected zero-momentum `tr F²` correlator gives `f₁²_bare = 7.14 ± 0.42` at seventeen sigma (bare level 0, `m_G = 1.51 ± 0.15 a`); every smearing level is strongly positive. This is [E], not [T]. At `β = 1.0` the signal is `0.2σ` and the residue is not extractable — because strong coupling is a coarse lattice where `m·a ≳ 3`, so the state is at or beyond the cutoff on `T = 6`. This is **not** evidence against P3: P3 is an analytic statement in the `β → 0` regime where `m·a` diverges, so the state is provably present but numerically sub-lattice-spacing. Strong-coupling provability and numerical positivity live at opposite ends of the coupling.
- **The `k → 0` curvature probe is positive but is not the `0⁺⁺` object.** On a single `L = 8` volume, the SU(3) `β = 5.7` softest nonzero mode `k = 2π/L = π/4` gives `G″ = +0.307 ± 0.055` in 6/6 transverse plane-wave modes, carried by `W_geo = +0.319 ± 0.005` (volume- and `k`-independent to 0.5%), with the action bracket at equipartition (`+0.037 ± 0.099`). SU(2) `β = 2.4` reproduces it (`+0.192 ± 0.058`, `W_geo = +0.252`). **But this is a `κ → 0` sign probe on a gauge-dependent gluon mode in a flat chart, not the gauge-invariant `0⁺⁺` residue.** Keeping only the proven `W_geo ≥ 0` and treating the rest statistically gives 95% one-sided lower bounds of `−0.098` (SU(3)) and `−0.035` (SU(2)); **zero is not excluded.** Closing it is one variance inequality, `(◆)` = H1 — not a better dataset. The gauge-invariant object is `m_GI = inf supp ρ_s` of `s = (β/2g) tr F²`.
- **The volume confound is removed.** At fixed `k = π/2`, `L = 4` and `L = 8` agree (`+0.425 ± 0.084` versus `+0.396 ± 0.063`); at fixed `L = 8`, `G″(π/4) = +0.307 < G″(π/2) = +0.396` but does not collapse. The earlier `L = 4 → 6` drop mixed a `k` decrease with a volume change; it is now separated. `W_geo` is pinned at `+0.32` across `L = 4, 6, 8` and both `k`.

### 9.4 The wall, located and instrumented

The honest image is not "we failed." It is that a shapeless wall was cut into a named object with a coordinate, a ruler, and an escape hatch. The wall is the **Savvidy concavity** (`Hess −D² + 2 ad F_B` indefinite; the negative transverse-magnetic mode survives the gauge quotient, measured 7.1% negative after `Q⊥`). The key is a **cutoff-consistent `G″(0) > 0` / coarse-mode gap lower bound** across `[1/48, 2.5]`, equivalently `f₁² > 0` transport, equivalently the gluon condensate `⟨¼ tr F²⟩` — **the same object under several names.** The residual is the one species: a one-sided variance upper bound at UV (L3) and IR (H1).

The wall did not fall. It became **transparent**: we can see the atom through the determinant, measure the sign through `W_geo`, and locate the failing cell exactly. It became **ruled**: roughly ninety routes, each with its precise death. It became **measured**: `f₁² > 0` at 17σ where resolvable, `D/S ≈ 2%`, `W_geo` constant to 0.5%, and a rigorous-style bound that is close to zero but not past it.

### 9.5 What would change the verdict

Only one thing: **the construction of the missing species** (or a proof that it cannot exist). Concretely, any one of:

- a nonperturbative proof of `Var_μ(⟨v,∇S⟩|0) ≤ E[H|0] + W − δ`, `(L, a)`-uniform (H1);
- a derivation of the deterministic negative-association inequality `Cov(G, Φ′|0) ≤ 0` for the frustrated SU(N) measure;
- a proof of the tachyonic-compensation inequality (TCI) pointwise and uniformly, hence convexity-at-the-vacuum of the nonperturbative effective potential plus absence of the Savvidy imaginary part;
- a nonperturbative transport of the strong-coupling gap/residue across `[1/48, 2.5]` with `a → 0`;
- a proof of the isolation I2 with a residue floor `f₁² ≥ r₀`, uniformly.

Each is a different name for the same missing machine. Until one exists, the Clay verdict is 0%, and it is honest to say so.

### 9.6 Sharpest line

> The wall did not fall; it became transparent, ruled, and instrumented. After roughly ninety independent frontal attacks, every honest route ends at the same sentence — an upper bound on a conditioned variance — written twice, once in the ultraviolet and once in the infrared, and the whole distance between this campaign and Clay is that one sentence. We proved the strong-coupling residue is positive, the lightest `0⁺⁺` cannot dissolve, the Haar curvature is non-negative, and the gauge-invariant vacuum is Perron–Frobenius; we proved a conditional theorem whose gluing lemmas are machine-checked; and we proved that positivity, convexity, spectral theory, Perron–Frobenius, and the renormalization group — the five exhausted paradigms — cannot supply the bound, because at the free point they all return zero. What remains is not a missing constant or a sharper estimate; it is a species that does not yet exist in mathematics: a single-sided, nonperturbative, non-circular, free-field-failing, cutoff-uniform upper bound on a variance. Build that, and the ultraviolet existence and the infrared gap close together. Until then, the honest word is not "proved" — it is "located."

---

## 10. References

Citations below reproduce the strings used in the corpus (author, venue, volume, year, page — no expanded titles are added). Where the corpus flags a bibliographic uncertainty, that flag is preserved rather than repaired.

1. A. Jaffe, E. Witten, *Quantum Yang–Mills Theory*, Clay Mathematics Institute. (Problem statement quoted in `docs/yang-mills-prior-work.md` §1.)
2. Streater–Wightman (1964), as cited in the Clay statement.
3. Osterwalder–Schrader, Comm. Math. Phys. **31** (1973) 83.
4. Osterwalder–Seiler, Ann. Phys. **110** (1978) 440.
5. Shen–Zhu–Zhu, Comm. Math. Phys. **400** (2023) 805, arXiv:2204.12737.
6. Kotecký–Preiss, Comm. Math. Phys. **103** (1986) 19.
7. Seiler, Lecture Notes in Physics **159** (1982).
8. Brydges–Giffen–Durhuus–Fröhlich, Nucl. Phys. B **275** (1986) 459.
9. Glimm–Jaffe–Spencer, Ann. Math. **100** (1974) 585–632.
10. Feldman–Osterwalder, Ann. Physics **97** (1976) 80–135.
11. T. Balaban, *Renormalization group approach to lattice gauge field theories* I–IV, Comm. Math. Phys. (1980s); volumes cited in the corpus as **109** (1987) 249; **116** (1988) 1; **122** (1989) 175, 355 (the corpus also records 102 (1985), 109 (1986), 122 (1989) and marks the volume numbers "to be verified" in one file).
12. Magnen–Rivasseau–Sénéor, Comm. Math. Phys. **155** (1993) 325–383 (one corpus file records the page range 417–430; the range is unverified here).
13. Källén, *Helv. Phys. Acta* **25** (1952) 417; Lehmann, *Nuovo Cim.* **11** (1954) 342.
14. Prékopa (1971).
15. Brascamp–Lieb, J. Funct. Anal. **22** (1976) 366.
16. Savvidy, Phys. Lett. B **71** (1977) 133.
17. Nielsen–Olesen, Nucl. Phys. B **144** (1978) 376.
18. Leutwyler, Nucl. Phys. B **179** (1981) 129.
19. Lüscher, Comm. Math. Phys. **54** (1977) 283.
20. Singer, Comm. Math. Phys. **60** (1978) 7.
21. Kato, *Perturbation Theory for Linear Operators* (1976), VII §3.
22. Reed–Simon IV, §XIII.1.
23. von Neumann–Wigner (1929).
24. Morningstar–Peardon, hep-lat/9901004.
25. Bauerschmidt–Brydges–Slade, Lecture Notes in Mathematics **2242** (2019).
26. Sinclair, Combin. Probab. Comput. **1** (1992) 351.
27. Jerrum–Sinclair, SIAM J. Comput. **18** (1989) 1149.
28. Feller, Ann. of Math. **55** (1952) 468.
29. Griffiths, J. Math. Phys. **8** (1967) 478.
30. Kelly–Sherman, J. Math. Phys. **9** (1968) 466.
31. Fortuin–Kasteleyn–Ginibre, Physica **57** (1972) 536.
32. Ginibre, Comm. Math. Phys. **16** (1970) 310.
33. Efron–Stein, Ann. Statist. **9** (1981) 586.
34. Royen, Far East J. Theor. Stat. **48** (2014) 139.
35. von Renesse–Sturm, Comm. Pure Appl. Math. **58** (2005) 923.
36. Guo–Shamai–Verdú, IEEE Trans. Inform. Theory **51** (2005) 1261.
37. Borcea–Brändén–Liggett, J. Amer. Math. Soc. **32** (2019) 205.
38. Cornwall, Phys. Rev. D **26** (1982) 1453.
39. Temple (1928); Temple–Weinstein / Temple–Kato as cited in the corpus.
40. Göpfert–Mack, Comm. Math. Phys. **81** (1981) 97.

---

### Appendix: verification map

| Claim in this paper | Label | File containing the proof / source |
|---|---|---|
| P3 residue theorem (I1) | [T] | `research/unconditional/PROOF-P3-UNCONDITIONAL.md`, `PROOF-P1P2.md`, `PROOF-P3.md` |
| P3 gives I1 not I2 | correction | `PROOF-P3-UNCONDITIONAL.md` §5; `PROOF-P1P2.md` §3 |
| No-dissolution theorem | [C] | `research/unconditional/PROOF-NO-DISSOLUTION-KINEMATIC.md`, `PROOF-DISSOLUTION-SPECTRAL-DENSITY.md` |
| Operative threshold `≈1.5 m_g` | [T] | `PROOF-ISOLATION.md` §4; `PROOF-EDGE-E2.md` |
| `W_geo ≥ 0`, `E[Φ″] ≥ N/12` | [T] | `research/unconditional/PROOF-JACOBIAN-SAVES.md` |
| Gauge-invariant PF | [T] finite V | `research/unconditional/PROOF-GAUSS-INVARIANT-PF.md` |
| Lemma SC | [T] | `research/unconditional/PROOF-COARSE-POINCARE.md` §(c) |
| Feshbach–Bogoliubov lemma | [T] abstract | `research/unconditional/FORGE-BLIND.md` §2; `D-BBS-INTEGRATE-NEG.md` (Lemma NS) |
| `G″ = S − D` | [T] identity | `INVENT-VARIANCE-INEQ.md`; `docs/new-tools-registry.md` A2 |
| Residue–level-gap lemma | [T] | `GAP-G2-edge`; `G2-ASSEMBLE-100.md` §2 |
| Conditional theorem `IN-0…IN-11` | [C] | `docs/conditional-theorem.md`, `lean/YMStatement/` |
| Spectral squeeze, H1 necessary-not-sufficient | [T] | `research/unconditional/REDUCTION-AUDIT.md` |
| Correct sufficient condition | [T] | `REDUCTION-AUDIT.md`; `LOWER-AMPLITUDE-MECHANISMS.md` §1 |
| Node board / honest % | [T]/[C]/[M] | `docs/PROOF-SKELETON.md`; `research/unconditional/G2-ASSEMBLE-100.md`, `G1-ASSEMBLE-100.md`, `ABC-FINAL-LEDGER.md` |
| Six tools and acid tests | [T]/[C] | `docs/new-tools-registry.md`; `INVENT-TRANSPORT-INEQ.md`, `INVENT-VARIANCE-INEQ.md`, `INVENT-NEW-IDENTITY.md`, `INVENT-GAUGE-INVARIANT-MASS.md`, `INVENT-NEW-MONOTONE.md` |
| No-gos | [N] | `INVENT-CONVEXIFICATION.md`, `INVENT-SECOND-CONE.md`, `INVENT-NEW-INDEX.md`, `SCHUR-CONVEX.md`, `GOLDS-POS-MASS.md`, `FORGE-FESHBACH-ANGLE.md`, `FORGE-KREIN-SPECTRAL.md` |
| One species (UV/IR), five conditions | [M] | `docs/PROOF-SKELETON.md` §3; `docs/YANG-MILLS-MASTER.md` §6 |
| Immune map | [N] | `docs/YANG-MILLS-MASTER.md` §5; `ABC-FINAL-LEDGER.md`; batch reports `ABC*`, `C-*`, `D-BBS-*`, `LOC-*`, `MS-*`, `LOWER-*`, `FORGE-*`, `K0-*`, `BS-*`, `TRANSMUTE-*` |
| `f₁² > 0` at β = 2.2 (17σ) | [E] | `research/unconditional/COMPUTE-F1SQ-SU2.md` |
| `G″` κ→0 probe (gauge-dependent mode) | [E] | `research/unconditional/COMPUTE-GPP-K0-CONTINUUM.md` |
| Dimensional correction `G″_a = a² m_phys²` | [T] | `docs/PROOF-SKELETON.md` (G1-CONTINUUM-SCALING) |
| "Machine that does not exist" | verdict | `docs/yang-mills-frontier.md` §7; `docs/YANG-MILLS-MASTER.md` §9; `TOOL-not-iii.md` |
