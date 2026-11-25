# T6-BINDING-ADVERSARIAL — separating **binding** from **isolation**: the program's discriminator is a **false dichotomy**; spectrally `I ⟹ B_spec ⟹ A` (binding is *implied by* isolation, not orthogonal to it), the `(binding, isolation)` square has an **empty quadrant**, and **binding + isolation still do not give the gap** (double well: gap→0, ρ→∞); the true independent axes are **atom-visibility `f₁²>0`** and **isolation `Λ★>m₁²`**, and the gap needs a **third, orthogonal input** (H1/ceiling + absolute scale)

> **2026-11-24 · Wave T6 — adversarial separation of binding and isolation.**
> **Target.** The R wave left a program: *binding is the discriminant of the atom, isolation is the gap*
> (`R5-EXACT-ANCHORS`, `R-WAVE-SYNTHESIS`), with R2's logical pair `I2 ⟹ I1`, `I1 ⇏ I2`, and R6's
> `H1 + fixed gap ⇏ atom`. The **next target** is a binding/isolation criterion at `β=2.5`. This note
> does the adversarial job: **separate the two, and try to kill the discriminator.** We construct the
> four cells (binding∧¬isolation, isolation∧¬binding, neither, both) on Källén–Lehmann measures, test
> them on double wells / Savvidy cells / 2D / 3D / strong coupling, and decide which is the true
> discriminant of **the atom** and which of **the gap**.
> **Reads:** `R-WAVE-SYNTHESIS.md`, `R2-ISOLATION.md`, `R5-EXACT-ANCHORS.md`, `R6-ATOM-ADVERSARIAL.md`,
> `Q6-H1-STATUS.md`; `docs/paper-v2/ch02-squeeze.tex` (squeeze, `eq:sufficient`, `eq:hstar`, `§schur`),
> `ch07-nogos.tex` (nested-cone, Schur, B2 single phase).
> **Witness:** `research/v3/w7/t6_binding_adversarial.py` (`py -3.12 -I`, numpy+scipy, deterministic,
> exit 0, **26/26 PASS**; double-well spectrum reproduces R6 §2.3 to printed digits).
> **Labels (charter §2).** **[T]** theorem/identity · **[C]** conditional · **[NG]** no-go for a named
> class · **[O]** exact missing step · **[V]** numerical witness. **No QED. Clay 0%.**

---

## 0 · Verdict (up front; no QED)

**The discriminator dies, and the split is the result.** Four findings:

1. **"Binding" is two inequivalent words, and the program conflates them `[T]`.**
   - **`B_spec` (spectral binding):** the bottom atom sits strictly **below the multiparticle
     continuum**, `f₁²>0` and `s_th := inf supp ν_ac > m₁²`. This is what "*stable single particle below
     the two-body threshold*" (R5) means on the spectral measure.
   - **`B_dyn` (dynamical binding):** a **specific mechanism** (the minimal connected `tr F²` polymer =
     one four-walled flux tube, P3) produces the pole. This is what "*the atom is a binding phenomenon*"
     (R5) means in the strong-coupling expansion.
   These are **not the same object**, and the program's "binding is the discriminant" is true of one and
   false of the other.

2. **Spectrally, binding is *implied by* isolation, so it is not an independent discriminant `[T]`.**
   Because `supp ν_ac ⊆ supp(ν − f₁²δ_{m₁²})`, one has `Λ★ ≤ s_th`, hence
   ```
        I (Λ★>m₁²)   ⟹   B_spec (s_th>m₁²)   ⟹   A (f₁²>0)               [T]
   ```
   **both strict.** So "binding" (spectral) is the *middle* of a sandwich: isolation implies it, it
   implies the atom. The task's requested cell **(b) "isolation but no binding" is EMPTY** under
   `B_spec` — there is no atomic measure that is isolated yet not a bound state below its own continuum.
   The `(binding, isolation)` square collapses to a **triangle**: `{¬B,¬I} ⊔ {B,¬I} ⊔ {B,I}`, with the
   `{¬B, I}` quadrant void. This is the direct kill of "binding vs isolation" as two independent levers.

3. **Even binding + isolation do not give the gap `[T]+[V]`.** The Savvidy double well keeps the atom
   `f₁²=O(1)` *and* a large isolation `λ₂−λ₁` while the transfer **gap Kramers-collapses** `λ₁~e^{−S}`
   and the mass ratio `ρ=√(λ₂/λ₁) → ∞` (witness; reproduces R6 §2.3). So `B ∧ I ⇏ G` — and the missing
   ingredient is exactly H1/ceiling, not more binding. The program's "isolation is the gap" is therefore
   `[C]` on **single-scale covariance with *bounded* `ρ`** (R2 §5), and the double well is the witness
   that the bounded-`ρ` side condition is **load-bearing, not benign.**

4. **The true discriminants are the two that are *not* binding.** The independent axes are
   **atom-visibility `A: f₁²>0`** (whose *mechanism* is `B_dyn` — R6's trace-anomaly/confinement
   physics) and **isolation `I: Λ★>m₁²`** (the second-eigenvalue/TP₂ datum, R2). **The gap `G: m₁²>0`
   is orthogonal to both** and needs a third input (H1/ceiling + an absolute scale). So the clean
   statement is:
   ```
        A  (⟸ B_dyn ;  mechanism = the bound state)          -- "which state is visible"
        I  (⟹ B_spec ⟹ A ;  stronger)                        -- "is the edge separated"
        G  (⊥ A, I without H1 + scale)                       -- "where is the floor"
   ```
   Binding is *not* a third axis: spectrally it is a function of `I`, dynamically it *is* `A`.

**One line.** *The (binding, isolation) dichotomy is false: spectrally `I⟹B⟹A`, so binding is redundant
given isolation (the isolation-without-binding cell is empty); dynamically binding is the atom's
mechanism; and binding + isolation together still miss the gap (double well). The true discriminants are
atom-visibility and isolation — **not** binding — and the gap is a third, orthogonal datum. No QED. Clay
0%.*

---

## 1 · Definitions: three predicates, and the two meanings of "binding"

Throughout `ν ≥ 0` is the Källén–Lehmann measure of the gauge-invariant `0⁺⁺` composite `O=tr F²`, on
`[0,∞)`, with `C(t)=∫e^{−√s\,t}dν(s)`, `f₁²=ν({m₁²})`, `m₁²=inf supp ν`.

| symbol | definition | meaning |
|---|---|---|
| **`A`** (atom / R2's `I1`) | `f₁² = ν({m₁²}) > 0` | the lightest state is **visible** to `O` |
| **`I`** (isolation / R2's `I2`) | `Λ★ := inf supp(ν − f₁²δ_{m₁²}) > m₁²` | the edge atom is **separated** from all heavier support |
| **`G`** (gap) | `m₁² = inf supp ν > 0` | the **floor** of support is positive |
| **`B_spec`** (spectral binding) | `f₁²>0` **and** `s_th := inf supp ν_ac > m₁²` | the atom is a bound state **strictly below the continuum** |
| **`B_dyn`** (dynamical binding) | a named mechanism (flux-tube polymer, P3) makes the pole | the **causal** story of how the atom arises |

`ν_ac` is the absolutely continuous part (the scattering continuum); `s_th` is its start. **`B_spec` is
the spectral restatement of R5's "stable single particle below the two-body threshold":** below the
continuum ⟺ below the multiparticle cut. **`B_dyn` is R5's "the minimal connected `tr F²` polymer is one
flux tube."** The R wave uses both without separating them; §2–§3 show they have opposite behaviour
relative to `I`.

Two elementary identities used below (all `[T]`, exact):
- `supp ν_ac ⊆ supp(ν − f₁²δ_{m₁²})` (a.c. mass is disjoint from a point atom), hence **`Λ★ ≤ s_th`**.
- `A ⇔ f₁²>0`; `I` and `B_spec` each *contain* `A` by definition.

---

## 2 · The logical lattice: `I ⟹ B_spec ⟹ A`, with both arrows strict

> **Theorem T6-1 (`I ⟹ B_spec ⟹ A`) `[T]`.**
> (i) `I ⟹ B_spec`: `Λ★>m₁²` and `s_th ≥ Λ★` give `s_th>m₁²`; with `I⟹A` (R2 Lemma A) this is `B_spec`.
> (ii) `B_spec ⟹ A`: by definition. Both implications strict (T6-2, T6-3). ∎

So the R-wave picture "atom is downstream of isolation, binding is the discriminant" is **half right and
half self-defeating**: the atom *is* downstream of isolation (R2's correction), but **binding is also
downstream of isolation** — it is *between* `I` and `A`, so it cannot be an *independent* discriminant
of either. The chain is a line, not a fork:
```
    I  ========►  B_spec  ========►  A                    (all [T], both strict)
    │             (s_th>m₁²)        (f₁²>0)
    └─ isolation is the parent; binding is its shadow; the atom is the grandchild
```

> **Theorem T6-2 (`A ⇏ B_spec`; the atom is not a binding) `[T]`.**
> There is `ν` with `f₁²>0` but `s_th = m₁²`: an **atom sitting exactly at its own continuum threshold**.
> *Proof / witness.* `ν = f₁²δ_{m₁²} + 1_{(m₁²,∞)}ρ(s)ds` (a.c. density `ρ>0` immediately above the atom).
> Then `f₁²>0` yet `s_th = m₁²`. **This is a "threshold state" / half-bound resonance.** ∎
> **Consequence:** the atom is **strictly weaker** than binding; `B_spec` is a genuine strengthening.
> The witness (cell (c)) is machine-checked: `A=T, B_spec=F, I=F`.

> **Theorem T6-3 (`B_spec ⇏ I`; binding is not an isolation) `[T]`.**
> There is `ν` with `B_spec` but `Λ★ = m₁²`: a bound state below the continuum, but with **discrete
> excitations accumulating at the bottom atom**.
> *Proof / witness.* `ν = f₁²δ_{m₁²} + Σ_{n≥1}w_n δ_{m₁²+1/n} + 1_{[s_th,∞)}ρ`, `s_th>m₁²`. The tower
> `m₁²+1/n → m₁²` gives `Λ★=m₁²` while `s_th>m₁²`. **Witness:** isolation margin `Λ★−m₁² = 1/N → 0` as
> `N→∞` (machine-checked `1.0e-4 → 1.0e-5`), i.e. **no uniform lower bound** on the isolation. ∎
> **Consequence:** the strong-coupling residue theorem (P3) is exactly a *dynamical* instance of this
> cell — it binds the flux tube (`B_dyn`, `f₁²≥r₀>0`) and **provably gives no `I₂`** (R5: amplitude
> grading cannot resolve `λ₂−λ₁`). **The program's one rigorous binding anchor is a
> binding-without-isolation model.** `[T]`.

**The `(B_spec, I)` square is a triangle (not a full 2×2).** By T6-1, `I ∧ ¬B_spec` is empty:

|  | `¬I` | `I` |
|---|---|---|
| **`¬B_spec`** | **(c)** threshold atom / gapped atomless — "neither" | **∅** (T6-1: `I⟹B_spec`) |
| **`B_spec`** | **(a)** tower accumulation / **strong-coupling P3** — "binding, no isolation" | **(d)** ordinary bound state / QCD glueball — "both" |

Under `B_dyn` (not `B_spec`) the `¬B,I` quadrant *can* be populated, but only by an atom produced by a
**non-binding** mechanism — a carrier/symmetry/trace-anomaly overlap rather than a bound state. **The
corpus has no controlled pure-YM model of that** (§4, `[O]`): in every exact solvable model the atom is
either zero/at-`m=0` (2D, free) or produced by binding (P3). So:

> **No-go T6-NG1 (the discriminator kill) `[T]/[NG]`.**
> *In the class of Källén–Lehmann measures, "binding" is not a primitive independent of "isolation":*
> `I ⟹ B_spec ⟹ A`, both strict, and the `(isolation ∧ ¬binding)` cell is empty. Consequently **no
> biconditional of the form "binding ⟺ atom", "binding ⟺ isolation", or "isolation ⟺ gap" holds**,
> and the program's "binding is the discriminant / isolation is the gap" cannot be a *pair of
> orthogonal discriminants*. The only independent spectral axes are **`A`** and **`I`**; "binding" is a
> function of them (implied by `I`); "gap" is orthogonal to both.

---

## 3 · The gap is orthogonal to binding and to isolation (the sharpest no-go)

The R-wave claim "isolation is the gap modulo IN-11" (`R2 §5`, `R6 §4`). We test both directions.

> **Theorem T6-4 (`A ∧ G ⇏ I`; a *gapped atom* is not isolated) `[T]`.**
> `ν = f₁²δ_{m₁²} + 1_{(m₁²,∞)}ρ`, `m₁²>0`: atom (`A`), gap (`G`), but `Λ★=m₁²` (`¬I`).
> *Witness:* `m₁=2, A=T, G=T, I=F`. **Extends R2's `I1⇏I2` to the gapped case:** even an atom *with a
> genuine positive gap* need not be isolated. Isolation is a strictly independent, stronger datum. ∎

> **Theorem T6-5 (`I ⇏ G`; isolation is not the gap) `[T]`.**
> (i) *Logical:* `ν = f₁²δ_0 + 1_{[3,∞)}ρ`: a **massless** atom at the bottom, isolated above
> (`Λ★=3>0=m₁²`, so `I`) but `G` fails (`m₁²=0`). *Witness:* `m₁=0, I=T, G=F`.
> (ii) *Physical (the decisive one):* the **Savvidy double well** has `B_spec ∧ I` with the transfer gap
> Kramers-collapsing `λ₁~e^{−S}→0` and `ρ=√(λ₂/λ₁)→∞` (witness §6·B; reproduces R6). So **uniform `I`
> does not give `G`** — the single-scale covariance hypothesis `Λ★=ρ²m₁²` with **bounded `ρ`** (R2 §5)
> is exactly what fails. ∎

> **Theorem T6-6 (`B ⇏ G`; binding is never the gap) `[T]`.**
> `B_spec` with `m₁²=0` (massless atom below a continuum gap) → `¬G` (T6-5(i)). `B_dyn` with the double
> well → residue `O(1)`, gap→0 (T6-5(ii)). **Neither the spectral nor the dynamical binding gives the
> gap, in either direction, without an independent H1/ceiling + absolute-scale input.** ∎

**The double-well witness (reproducing R6 §2.3).** From the witness (§6·B), transfer/SUSY operator
`H=−∂²+W`, `W=V'²/4−V''/2`, `V=(g/2)(x²−1)²`:

```
   g      λ1=GAP      f₁²       I₂=λ2−λ1     ρ=√(λ2/λ1)
  0.5   7.921e-01   1.0144       2.76          2.12
  2.0   7.487e-01   0.8010       5.25          2.83
  8.0   1.191e-01   0.8842      13.12         10.54
 12.0   2.484e-02   0.9249      20.01         28.39
```
The residue and the isolation **grow** while the gap **collapses** and `ρ` **diverges**. This is the
exact refutation of "binding + isolation ⟹ gap" (uniformly): the missing datum is H1, not binding. It is
also R6's "residue floor is the tractable half" made structural — the **binding + isolation are both
robust, the gap is the fine-tuned object.**

---

## 4 · Model battery: double wells, Savvidy, 2D, 3D, strong coupling

| model | `A` | `B_spec` | `B_dyn` | `I` | `G` | cell | label |
|---|---|---|---|---|---|---|---|
| **strong coupling P3** (`β<1/48`) | ✓ `f₁²≥r₀` | ✓ (below cut) | **✓ flux tube** | **✗** (R5: no `I₂`) | separate (OS/SZZ) | **(a)** | **[T]** |
| **Savvidy double well** | ✓ `O(1)` | ✓ (no continuum) | ✓ | ✓ large | **✗** `λ₁→0`, `ρ→∞` | **(d) w/o gap** | **[T]+[V]** |
| **2D YM, `tr F²` carrier** | atom at `m=0` | ✗ | ✗ (`tr F²=H` conserved) | ✗ | ✗ | (c)/degenerate | **[T]** |
| **2D YM, holonomy carrier** | — | — | — | ✓ `ρ²=8/3` | ✓ | (different `ν`!) | **[T]** |
| **3D / 4D lattice** | ✓ `[V]` | ✓ | ✓ (confining) | ✓ `[V]` | ✓ `[V]` | (d) | **[T-standard]/[V]** |
| **free field / `:φ²:`** | ✗ `f₁²=0` | ✗ | ✗ | ✗ | ✓ (massive) | (c) | **[T]+[V]** |

**The 2D carrier mismatch — a correction to the program's evidence base `[T]`.** R2 §4 W2 quotes the 2D
isolation `ρ²=C₂(1)/C₂(½)=8/3` as "a symmetry-protected, exact positive isolation." But R5 §1.1 proves
that in 2D `tr F² = tr E² = H` is the **conserved energy density**, so `⟨0|tr F²|G⟩ = C₂(R)⟨χ_0|χ_R⟩ = 0`
— the `tr F²` atom is **zero** and the `tr F²` two-point function is a `t`-constant contact. The
`ρ²=8/3` tower is the **holonomy/character** spectrum, a *different operator*. Hence **R2's positive
isolation anchor and R5's zero atom anchor are statements about two different measures `ν`** — the
isolation R2 computes is not the isolation of the `tr F²` measure that `H*` (`ch02 eq:hstar`) actually
uses. Witness §6·D confirms `<χ_0|χ_R>=δ_{0R}` (to `1e-16`) and that `8/3` is the character ratio.
**Consequence:** the corpus's only *exact* positive isolation is for the wrong carrier; the `tr F²`
measure in 2D is degenerate (no atom at `m>0`, no isolation). This removes one of the R-wave's two
"exact isolation anchors" (`R-WAVE §1`) from the `tr F²` program. **[T]**

**Strong coupling is the corpus's own cell (a).** P3 (`PROOF-P3-UNCONDITIONAL`) proves `f₁²≥r₀(β)>0`
(`B_dyn`: the flux tube is a single bound particle) but **explicitly no `I₂`** (R5 §1.4: "the cost
grading grades polymer *amplitudes*, not transfer-matrix *eigenvalue spacing*"). The witness §6·C makes
this quantitative: a `t`-uniform `O(u)` remainder cannot resolve the isolation `f₂²e^{−(λ₂−λ₁)t}` until
`t ~ (1/(λ₂−λ₁))log(f₂²/(Bu))`, which diverges as the isolation shrinks — exactly the `Cu` vs
`Be^{−δt}` obstruction of `SUB2-I2-STRONG`. **So the program's binding theorem proves binding, not
isolation; the two are empirically (in the one rigorous case) in *different cells*.** `[T]`

---

## 5 · The strongest no-go, and the weakest hypotheses

### 5.1 The headline no-go (restated with its class)

> **T6-NG1 `[T]/[NG]` (class = Källén–Lehmann measures; and frameworks preserving RP + a confined
> single phase + the `tr F²` carrier).** Binding is not an independent discriminant:
> `I ⟹ B_spec ⟹ A` (both strict), the `(I ∧ ¬B_spec)` cell is empty, and **`B ∧ I ⇏ G` even with a
> genuine gap elsewhere** (double well). Therefore the program **cannot be organized as a pair
> (binding-discriminant, isolation-gap)**: binding is redundant given isolation (spectral) or *is* the
> atom (dynamical), and the gap needs the orthogonal H1/ceiling + scale. **Any "binding theorem" is,
> at best, a proof of `A` (or `B_spec`), never a proof of `I` and never a proof of `G`.**

The single sharpest *model-level* statement: **the strong-coupling residue theorem — the corpus's one
unconditional binding anchor — is a `B_dyn`-without-`I` object (cell (a)), and the Savvidy double well
is a `B∧I`-without-`G` object (cell (d) minus gap).** The two together exhaust the R-wave's positive
ledger and show it does not reach the gap by *either* route.

### 5.2 Weakest hypotheses (the task's explicit ask)

> **(H1) `B_spec ⟹ I`** needs exactly one hypothesis `[T]`:
> **(H-fin)** the discrete spectrum of `ν − f₁²δ_{m₁²}` **below `s_th` is finite** (no accumulation at
> the bottom atom). Then the finitely many levels `m₁²<m₂²<…<m_K²<s_th` give `Λ★=m₂²>m₁²`, i.e. `I`.
> *This is the corpus's finite-tower theorem (R2 W5 / Prony) in isolation form* — physically generic
> (a glueball tower below `2m_g` is finite), but an **extra hypothesis, not a consequence of binding.**

> **(H2) `I ⟹ B_spec`** needs **no hypothesis** `[T]` — it is T6-1, automatic. **This is the sharpest
> asymmetry: one direction of "binding ↔ isolation" is free, the other needs a finiteness input.**
> Hence "isolation but no binding" cannot exist; "binding but no isolation" needs an accumulating tower.

> **(H3) `B_dyn ⟹ I`** needs an **independent order-2 / total-positivity hypothesis** `[O]`, which is
> **not weaker than the isolation theorem itself** (`(E-I2)`, `ρ≥1+c`, R2 §6.1). P3 is the witness that
> `B_dyn` holds while `I` fails. So a *dynamical* binding theorem does **not** imply isolation short of
> assuming the isolation content.

> **(H4) `I ⟹ G`** needs **single-scale covariance with *bounded* `ρ`** `[C]` (`Λ★=ρ²m₁²`,
> `1<ρ≤ρ_max<∞`, R2 §5). The double well shows that **unbounded `ρ`** (`ρ→∞`, gap→0) breaks it; so the
> `[C]` is load-bearing. Without SSC, `I ⇏ G` (massless isolated atom).

> **(H5) `B ⟹ G`** has **no non-circular hypothesis** `[NG]`. Spectral binding holds with `m₁²=0`;
> dynamical binding holds in the double well with gap→0. To get `G` from binding one must assume
> `m₁²>0` — i.e. `G` itself. **Binding never gives the gap.**

**Reading.** The "weakest hypotheses" table is itself the kill: `H2` is free (isolation ⟹ binding), `H1`
needs a finiteness input, and `H3–H5` show binding never reaches the gap. The pair (binding, isolation)
is **linearly dependent** (`H2`), not a basis — the program's target as a *pair* is wrong; the target is
the single ordered pair **(`A` mechanism, `I`)** plus the orthogonal ceiling **(H1-condition `χ≤1/λ`)**.

---

## 6 · Witness — `research/v3/w7/t6_binding_adversarial.py` (`py -3.12 -I`, exit 0)

Deterministic, numpy 2.4.6 + scipy 1.18.1, no MC, **26/26 PASS**, ~2 s. Quoted output (abridged):

```
PART 1  four cells + implication lattice
(a) tower: m1=1.0 f1=0.5 s_th=4.0 ; Lamb-m1: 1.0e-4 (N=1e4) -> 1.0e-5 (N=1e5)
    A=True B_spec=True G=True  ; isolation margin -> 0  (NO uniform I)
(b) forcing ~B_spec (s_th=m1) forces NOT I  -> cell (b) EMPTY under B_spec
(c) threshold atom: m1=1.0 f1=0.4 s_th=1.0 ; A=True B_spec=False I=False
    (c') gapped atomless: A=F B=F I=F G=T
(d) both: m1=1.0 f1=0.5 Lamb=2.5 s_th=4.0 ; A & B_spec & I & G all TRUE
random battery (4000 measures): I=>B_spec violations=0 ; B_spec=>A violations=0 ; I=>A violations=0
A&G but NOT I : m1=2.0 A=T G=T I=F        (T6-4)
I but NOT G   : m1=0.0 Lamb=3.0 I=T G=F   (T6-5)

PART 2  Savvidy / double well (reproduces R6 section 2.3)
   g      lam1=GAP     f1^2       I2=lam2-lam1   rho
  0.5   7.921e-01   1.0144       2.757          2.12
  2.0   7.487e-01   0.8010       5.255          2.83
  8.0   1.191e-01   0.8842      13.123         10.54
 12.0   2.484e-02   0.9249      20.006         28.39
  => residue O(1) and I2 large while GAP collapses and rho -> inf  (B & I not enough for G)

PART 3  strong coupling: t-uniform O(u) hides isolation until T ~ (1/(c m1)) log(f2^2/(Bu))
  c=0.02 -> T ~ 150 ; amplitude plateau reached ~ 1.4  (P3 grades amplitudes, not lam2-lam1)

PART 4  2D carrier mismatch
  <chi_0|chi_R> = 1.000000, 2.6e-16, -2.7e-16   (orthogonality => <0|trF^2|G>=0)
  C2(1)/C2(1/2) = 8/3  = the HOLONOMY ratio, not a tr F^2 atom property
```

Readings: **(a)** binding with isolation margin `→0` (no uniform `I`; strong-coupling pattern);
**(b)** forced-empty — isolation-witout-binding does not exist spectrally; **(c)** atom with
`B_spec=F`, `I=F` (threshold state); **(d)** both; **double well** residue+isolation robust while the
gap collapses; **2D** the `8/3` isolation is the holonomy's, the `tr F²` atom is zero. Witness, not
proof; the theorems are §2–§5, the numerics are the structural illustrations.

---

## 7 · Falsification battery for any proposed **binding theorem** (FB1–FB6)

A proposed binding theorem proves, uniformly in `(L,a)` at `β=2.5`, `f₁²≥r₀>0` (and perhaps claims `I₂`
or `G`). It must survive **all six**; failing any refutes it. Each is an exact/machine-checked
countermodel in the corpus or this note.

- **FB1 — threshold / BIC (must distinguish *below-continuum* from *embedded*) `[T]+[V]`.**
  Feed cell (c) `ν=f₁²δ_{m₁²}+1_{(m₁²,∞)}ρ` (atom at `s_th`). A "binding" theorem that fires here is
  really an *atom* theorem — it cannot tell a stable bound state from an embedded/threshold state.
  *Refutes any theorem that conflates `A` with `B_spec`.*
- **FB2 — tower accumulation (must NOT infer isolation from binding) `[T]+[V]`.**
  Feed cell (a) `ν=f₁²δ_{m₁²}+Σ_n w_nδ_{m₁²+1/n}` (`B_spec` true, `Λ★=m₁²`). A binding theorem claiming
  `I₂` is refuted — exactly P3's own boundary (`Cu` vs `Be^{−δt}`). *Refutes "binding ⟹ isolation".*
- **FB3 — double well / Savvidy (must return gap→0 and ρ→∞) `[T]+[V]`.**
  Feed `V_eff=(g/2)(x²−1)²`: `f₁²=O(1)`, `I₂` large, `λ₁~e^{−S}→0`, `ρ→∞`. A binding theorem claiming
  the **gap** or a **uniform ρ** is refuted. *Refutes "binding ⟹ gap" and any uniform-ρ isolation.*
- **FB4 — free field / `:φ²:` (must return `f₁²=0`) `[T]+[V]`.**
  A theorem firing on the free massive composite (`ν` = two-particle cut, `f₁²=0`) is **vacuous**
  (Schur-concave obstruction, `ch02 §schur`). *Acid test.*
- **FB5 — carrier + channel specificity (must be operator-labelled) `[T]+[V]`.**
  On the 2D **holonomy** carrier the solution has a discrete `ρ²=8/3` tower; on the **`tr F²`** carrier
  the same 2D theory has `f₁²=0` (`tr F²=H` conserved). A binding theorem that is **carrier-blind**
  wrongly assigns the holonomy's binding/isolation to `tr F²` (the §4 mismatch). *Refutes any theorem
  stated for "the `0⁺⁺` channel" without naming `O` (R6 carrier hazard).*
- **FB6 — `Z_O` / dimension (floor must be a *ratio*) `[T]`.**
  Under `f₁²(μ)=Z_O²f₁²(μ₀)` (`Z_O²>0` a square) `0→0` but `r₀` rescales; and `G″_a=a²m²`. A theorem
  producing an **absolute** `r₀` without `Z_O²`-covariance or `a²`-scaling is scheme/dimension-wrong
  (R6-c; Q6-F5). *The honest floor is `f₁²/χ≥c` or `f₁²/M₀≥c`, a ratio.*

**What a valid binding theorem therefore is.** A **carrier- and channel-specific**, `Z_O²`-covariant,
`a²`-scaled proof of `f₁²≥r₀` that **returns `0` on the free composite** (FB4), **does not fire on
threshold states** (FB1), **does not claim isolation** (FB2) **or the gap** (FB3) — i.e. it proves
**`B_dyn` / `A`**, exactly what P3 does at `β<1/48`. It is a genuine sub-target but **not** the gap and
**not** isolation.

---

## 8 · Status ledger

| item | statement | label |
|---|---|---|
| spectral binding `B_spec` | `f₁²>0 ∧ s_th>m₁²` (atom below a.c. continuum) | def. |
| dynamical binding `B_dyn` | flux-tube polymer makes the pole (P3, R5) | def. |
| **T6-1** `I ⟹ B_spec ⟹ A` (both strict) | `Λ★ ≤ s_th`; sandwich | **[T]** |
| **T6-2** `A ⇏ B_spec` | threshold atom `s_th=m₁²`, `f₁²>0` | **[T]+[V]** |
| **T6-3** `B_spec ⇏ I` | tower `Λ★=m₁²` below `s_th`; margin `1/N→0` | **[T]+[V]** |
| **T6-4** `A ∧ G ⇏ I` | gapped atom touching continuum (extends R2's `I1⇏I2`) | **[T]+[V]** |
| **T6-5** `I ⇏ G` | massless isolated atom; **double well `ρ→∞`, gap→0** | **[T]+[V]** |
| **T6-6** `B ⇏ G` (either sense) | `B_spec` massless; `B_dyn` double well | **[T]+[V]** |
| **T6-NG1** binding not an independent discriminant | `(I∧¬B_spec)` empty; `B∧I⇏G` | **[T]/[NG]** |
| cell (a) binding, no isolation | tower; **strong-coupling P3** | **[T]** |
| cell (b) isolation, no binding (`B_spec`) | **∅ — impossible** | **[T]** |
| cell (c) neither | threshold atom; free/`:φ²:` | **[T]+[V]** |
| cell (d) both | ordinary bound state; QCD; double well (w/o gap) | **[T-standard]+[V]** |
| 2D `ρ²=8/3` isolation | **holonomy** ratio, **not** the `tr F²` atom (which is `0`) | **[T]** (correction) |
| weakest hyp. `B_spec⟹I` | (H-fin) finite discrete tower below `s_th` | **[T]** |
| weakest hyp. `I⟹B_spec` | **none** (automatic) | **[T]** |
| weakest hyp. `B_dyn⟹I` | order-2/TP₂ `(E-I2)` — not weaker than isolation | **[O]** |
| weakest hyp. `I⟹G` | single-scale, bounded `ρ` (R2 §5) | **[C]** |
| weakest hyp. `B⟹G` | none non-circular | **[NG]** |
| true atom discriminant | `A: f₁²>0`, mechanism `B_dyn`, **not** binding-as-axis | — |
| true gap discriminant | H1/ceiling + absolute scale; **⊥ binding, ⊥ isolation** | `H*` = gap |
| free-field acid | `f₁²=0` for `tr F²`, `:φ²:` | **[T]+[V]** |
| **Clay mass gap** | unbuilt | **[O] · 0%** |

---

## 9 · Sharpest line

**The program's discriminator is a false dichotomy.** The R wave split "binding (the atom's discriminant)"
from "isolation (the gap)" as if they were independent levers. They are not. On the Källén–Lehmann
measure the two are **nested**: `supp ν_ac ⊆ supp(ν−f₁²δ)`, hence `Λ★ ≤ s_th`, hence
**`I ⟹ B_spec ⟹ A`**, both strict. So *spectral binding is implied by isolation* — the
"isolation-but-no-binding" cell is **empty** (a theorem, witness-confirmed), and binding cannot be an
independent discriminant of the atom; only the **dynamical** binding (the flux-tube mechanism, P3) is a
valid sufficient story for the atom, and it is exactly the atom itself (`f₁²>0`), giving **no isolation**
(the corpus's own cell (a)). In the other direction, **neither binding nor isolation gives the gap**:
the double well keeps `f₁²=O(1)` and a large `λ₂−λ₁` while the gap Kramers-collapses and `ρ→∞` — so
"binding + isolation ⟹ gap" is false, and "isolation is the gap" is `[C]` on single-scale covariance
with **bounded `ρ`** (a load-bearing, not benign, side condition). The two *truly independent* spectral
axes are **atom-visibility `f₁²>0`** (mechanism: binding/trace-anomaly overlap, R6) and **isolation
`Λ★>m₁²`** (the order-2/TP₂ datum, R2); the gap is a **third, orthogonal** datum (H1/ceiling + absolute
scale), which is why the full `H*` is the gap. **The target is not a binding/isolation pair; it is the
ordered pair (`f₁²>0`, `Λ★>m₁²`) plus the ceiling `χ≤1/λ`.** **No QED. Clay 0%.**

---

*Refs.* **R wave:** `R-WAVE-SYNTHESIS.md` (§1 isolation is parent, §2 binding-discriminant, §4
adversarial); `R2-ISOLATION.md` (Lemma A `I2⟹I1`, Lemma B `I1⇏I2`, Lemma C/D, SSC §5, `(E-I2)` §6.1,
W1–W6); `R5-EXACT-ANCHORS.md` (§1.1 2D `trF²=H` ⇒ atom 0; §1.4 P3 no `I₂`; §4 `(T★)`); `R6-ATOM-ADVERSARIAL.md`
(§2.3 double well residue robust/gap collapses; §4 R6-a gapped atomless; §4 R6-c `Z_O`). **Paper:**
`docs/paper-v2/ch02-squeeze.tex` (squeeze `eq:squeeze`, thin atom, `§schur`, `eq:sufficient`,
`eq:hstar`); `ch07-nogos.tex` (Schur master, nested-cone, B2 single phase). **Q6:** `Q6-H1-STATUS.md`
(H1 strictly weaker; thin atom; double well; F1–F5). **Strong coupling:** `unconditional/PROOF-P3-UNCONDITIONAL.md`
(flux-tube residue, `β<1/48`, no `I₂`); `SUB2-I2-STRONG.md` (`Cu` vs `Be^{−δt}`). **Isolation banked:**
`G2D-ISOLATION-I2.md` (SSC, `ρ`), `PROOF-EDGE-E2.md` (order-2 blindness). **Standard (textbook, used as
such):** Källén (1952)/Lehmann (1954); Osterwalder–Schrader (1973–75); Migdal (1975) / Witten CMP 141
(1991) (2D YM exact solvability); Brascamp–Lieb (1976); Kato (1976) Ch. VII; trace anomaly
`Θ^μ_μ=(β/2g)tr F²`; Weyl integration formula (SU(2) characters). **Witness:**
`research/v3/w7/t6_binding_adversarial.py` (`py -3.12 -I`, numpy+scipy, deterministic, 26/26 PASS,
double-well reproduces R6 §2.3). **No third-party downloads; no `[C]` upgraded to `[T]`; no QED; Clay 0%.**
