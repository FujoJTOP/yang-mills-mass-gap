# R3-ATOM-FUNCTIONAL — building a genuinely discontinuous atom-reading functional: **many exist and are exact, but every one collapses onto the same weak-\* no-go through one of three equivalent faces (location / limit / family-assumption); none is computable from finite `(L,a)`-data without presupposing the gap**

> **2026-11-24 · Wave R3 — the atom, the irreducible target left by H1.**
> Target: `f₁² = |⟨0|tr F²|G⟩|² = ν({m₁²}) ≥ r₀ > 0`, `m₁² = inf supp ν`.
> The corpus's claim (`docs/paper-v2/ch09-weakstar.tex`, `FALSIFY-TOPOLOGICAL.md`): the atom is
> **exactly** the weak-\* discontinuity of the boundary density; every finite datum is weak-\* continuous;
> atomless measures are weak-\* dense; hence **no continuous functional forces the atom**
> (Thm 4, ch09). The task: build a genuinely **discontinuous** atom-reading functional out of the
> standard analytic machines (spectral projectors, Hankel/Toeplitz, de Branges/spectral shift,
> resolvent boundary values, Stone–Weierstrass), decide whether any is computable from the
> physical measure **without already knowing the gap**, and whether it evades the no-go or collapses.
>
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional/assumption · **[NG]** no-go
> for a named class · **[O]** exact missing step · **[V]** numerical witness. **No QED. Clay 0%.**
> **Witness:** `research/v3/w7/r3_atom_functional.py` (`py -3.12 -I`, numpy/scipy only, one job, no downloads).

---

## 0 · Verdict (up front; no QED, Clay 0%)

1. **Discontinuous atom-reading functionals EXIST and are EXACT `[T]`.** Four independent ones were
   built and machine-tested: the Stieltjes boundary jump `[C1]`, the sup-over-location atom
   functional `[C2]`, the moment-Hankel/Christoffel reading `[C5]`, and the finite-moment blindness
   witness `[C3]`. Each returns `f₁²` exactly and passes the acid test (returns `0` on atomless
   measures). So the answer to "does a discontinuous atom-reading functional exist?" is **yes, trivially,
   many** — the atom functional itself is one.

2. **But every one is circular, in one of exactly THREE equivalent faces `[T]`/`[NG]`:**
   - **(R1) location** — the functional must be evaluated *at* `x = m₁²`; setting `x` presupposes the gap.
   - **(R2) limit** — evaluating it needs `d,ε → 0` or `n → ∞` (the full measure, **infinite data** = A3).
   - **(R3) family** — the only location-free variant `[C2]` targets the *wrong* atom (the largest, not
     the edge), and is positive on any **finite-volume** free field.
   **R1 = R2 = R3 are one obstruction**, three coordinates on the weak-\* discontinuity: atom = a `δ` in
   the boundary density; finite data = bounded-kernel integrals = continuous; atomless = dense.

3. **Answer to "computable without already knowing the gap?": YES given the exact measure, NO given any
   finite resolution — and finite resolution is all physics supplies `[T]`/`[NG]`.** From the full `ν`
   the edge atom is computable by the double limit `lim_{d↓0}lim_{ε↓0}` with **no knowledge of `m₁²`**
   (you look at the lowest spike). But a **finite-data** instrument `g` cannot resolve it: for any
   weak-\* continuous `g`, `sup_ν |E(ν) − g(ν)| ≥ f₁²/2` (**Stone–Weierstrass resolution bound `[T]`**),
   i.e. the **instrument's error is at least the signal** — circular (Theorem R3-2, §3).

4. **The exact missing input is unchanged and unmoved `[O]`:** a uniform-in-`(L,a)` **resolution** bound
   at the edge — resolve below the isolation gap `I₂` — which is the atom floor *with isolation*, i.e.
   `H*` itself, not a consequence of any functional. A discontinuous atom-reader is a **reading**, and
   its uniform evaluation **is** the target.

5. **One genuine side-finding (a correction) `[T]`:** the literal ch09 atom-weight formula with the
   window set equal to the kernel width (`d = ε`) returns `f₁²/2`, not `f₁²` — a factor-2 slip. The
   correct object is the **double limit** (§4, `[C1]` table; `[T]`).

---

## 1 · The target, restated as an edge-atom functional

The Källén–Lehmann measure `ν ≥ 0` on `[0,∞)` of `tr F²` has bottom `m₁² = inf supp ν` and edge atom

```
f₁² = E(ν) := ν({ m₁² }) = ν({ inf supp ν }) .                                            (1.1)
```

**The composition observation `[T]`.** Write `E = ev ∘ b` where `b(ν) = inf supp ν` (the bottom) and
`ev_x(ν) = ν({x})` (the atom mass at `x`). **Both `b` and `ev_x` are weak-\* discontinuous, and on the
relevant sequence they fail together:**

- `b`: `ν_h = (1−w)δ_{s₂} + w·Unif[m−h, m+h] → ν = (1−w)δ_{s₂} + wδ_m` weak-\*, but
  `b(ν_h) = m − h → m = b(ν)` **from below** (one-sided).
- `ev`: along the same sequence, `ev_m(ν_h) = 0` but `ev_m(ν) = w` — a genuine **jump of size `w`**.

So `E` is *doubly* discontinuous. This is not a pathology but the definition of the atom: the edge
atom is the boundary value of the density **at a point the measure has to be resolved to locate**. The
rest of this file shows every analytic machine that can express `E` reproduces this, and that the
discontinuity cannot be "cashed in" at finite resolution.

---

## 2 · The five named candidate classes — built, tested, collapsed

We treat each machine of the task, give its **exact** form for `E`, its **acid test**, and its
**collapse** onto R1/R2/R3. Witness numbers in §4.

### (a) Spectral projectors from OS reconstruction — collapses to (R1) or to rigid multiplicity `[T]`/`[NG]`

In OS reconstruction `f₁² = ⟨0| O P_{Δ} O |0⟩` where `P_Δ = 1_{Δ}(H)` is a spectral projector of the
physical Hamiltonian and `Δ ∋ m₁²` isolates the edge atom. Two variants, both collapse:

- **Interval projector `P_{[0,Λ]}`.** `⟨0|O P_{[0,Λ]} O|0⟩ = ν([0,Λ]) → E(ν)` only as `Λ ↓ m₁²`. To
  isolate the atom you must **choose `Λ = m₁² + I₂`** — you must know the gap `m₁²` and the isolation
  `I₂`. That is **(R1)**: the projector needs the eigenvalue it is supposed to detect.
- **Point projector `P_{\{m₁²\}}`** (rank-1 onto the edge eigenvector). Its **dimension** is the
  multiplicity — an **integer**, hence rigid (the corpus's no-go B3 / `FALSIFY-FLOW`: a signed index is
  an integer, equal to its free value). The **residue weight** `f₁²` is an **overlap** `|⟨G|O|0⟩|²`, not
  a count, so the projector's rank is **structurally blind to `f₁²`**. `[NG]`.

The positive content `[T]`: `E(ν)` *is* an interval-projector expectation in the limit `Λ↓m₁²`; the OS
machine reproduces `[C1]` exactly. It does not evade R1.

### (b) Hankel / Toeplitz — collapses to (R1)+(R2) `[T]`/`[NG]`

The moment sequence `M_k = ∫s^k dν` is **injective** (Hausdorff/moment problem) — so `E` **is** a
function of the full `(M_k)` (ch09 §topology says this). The Hankel matrix `G_n = (M_{j+k})` and the
reproducing kernel give the **Christoffel function**

```
λ_n(x) = 1 / (vₙ(x)ᵀ G_n⁻¹ vₙ(x)),    vₙ(x) = (1, x, …, xⁿ)ᵀ .                            (2.1)
```

The classical asymptotics `[std]` (Geronimus): `liminf_n λ_n(x) ≥ ν({x})` for every `x`, and
`λ_n(x) → dν_ac/dx` at a.c. points. Hence at the **atom** `λ_n(m₁²) ↛ 0` (stays `≥ f₁²`; empirically
`→ f₁²`, witness `0.3001`), at the **gap** `λ_n → 0` (`0.005`), and in the **continuum** `λ_n → ` density.
So the Christoffel sequence reads `f₁²` — **but only** (i) **evaluated at `x = m₁²`** (**R1**) and
(ii) **as `n → ∞`** (**R2** = infinite data). No finite `n` certifies the atom (witness: `n = 4` gives
the atom `0.303` and the *gap* `0.246` — a gap point masquerades as an atom at finite `n`). `[NG]`.

Equivalent Toeplitz/Szegő and Schur/Pick readings (Christoffel `↑ f₁²`, the Nevanlinna–Pick `[0,f_max]`
range, Toeplitz determinants) all bound `f₁²` **from above** and need the same location/limit — the ch06
`TOOL-TOEPLITZ` / `TOOL-NEVANLINNA` rows, now seen as one object. `[T]`/`[NG]`.

### (c) de Branges / spectral shift — collapses to (rigid integer) or to (A3) `[T]`/`[NG]`

The **spectral shift function** `ξ(λ)` of a pair `(H, H₀)` (Birman–Krein) is a **step function**: it
**jumps by the eigenvalue multiplicity** at a bound state, flat below the bottom of the essential
spectrum. So `ξ(λ)` is a genuine discontinuous atom-**presence** detector. But:

- `ξ` carries the **multiplicity** (an integer) — **rigid** (B3), blind to the residue `f₁²` (the
  overlap), exactly as in (a). `[NG]` for the weight.
- The **weight** lives in the de Branges **characteristic function** `E(z)` (the transfer operator's
  `det`): an edge bound state is a **Blaschke factor** of `E` at the eigenvalue. Reading off the
  Blaschke mass = the residue requires the **Akhiezer/scattering factorization** at the (unknown) real
  point `m₁²` — i.e. **boundary values on the real axis** = **(R1)+(R2) = A3**.

So the de Branges machine reproduces R1+R2 (or the rigid-integer no-go). `[NG]`.

### (d) Resolvent boundary values — the canonical (R1)+(R2) `[T]`/`[NG]`, with a factor-2 correction

The **Stieltjes transform** `Π(z) = ∫ dν(s)/(s−z)` is Herglotz; `(1/π) Im Π(x+iε) = ∫ P_ε(x−s)dν(s)`
(`P_ε` = Poisson kernel, an approximate identity) `→ ν` weak-\* as `ε↓0`. The **exact** atom read is the
**double limit**

```
ν({x}) = lim_{d↓0} lim_{ε↓0} (1/π) ∫_{x−d}^{x+d} Im Π(u+iε) du .                           (2.2)
```

**Correction `[T]` (charter §2).** `ch09` eq. `atomweight` prints (2.2) with the window set to `d = ε`.
With `d = ε` the inner integral of the Poisson mass over `[x−ε, x+ε]` is `(1/π)[arctan(1) − arctan(−1)]
= 1/2`, so the literal formula returns `f₁²/2`, **not** `f₁²`. The correct statement is the double limit
`d → 0` **after** `ε → 0` (witness §4 `[C1]`: atom → `0.300`, the `d=ε` variant → `0.150`). This is a
factor-2 slip of normalization, now recorded and fixed.

Witness `[C1]`: at the atom `x = m₁² = 1`, `M(d) → 0.30` for all `d`; at an atomless point `→ 0`. The read
is **exact**, but needs `x = m₁²` (**R1**) and `d, ε → 0` (**R2** = A3). **Collapses onto the no-go.**

### (e) Stone–Weierstrass — a quantitative resolution bound, closing model-free forcing `[T]`/`[NG]`

Stone–Weierstrass says continuous data (bounded-kernel integrals = finite moments/evaluations) is
**dense in `C(𝓜)`** — the *continuous* dual. A discontinuous `E` is **not** a uniform limit of continuous
data. Quantitatively:

> **Theorem R3-2 (resolution bound) `[T]`.** Let `g: 𝓜(K) → ℝ` be weak-\* continuous (any finite-data
> functional). For any `ν₀` with `E(ν₀) = w > 0`,
> ```
>   sup_{ν ∈ 𝓜(K)} |E(ν) − g(ν)|  ≥  w/2 .
> ```
> *Proof.* Let `ν_h` be the edge bump of `ν₀` (`E(ν_h)=0`, `ν_h → ν₀` weak-\*). Then `g(ν_h) → g(ν₀)` by
> continuity, and
> `|E(ν₀)−g(ν₀)| + |E(ν_h)−g(ν_h)| ≥ |E(ν₀)−E(ν_h)| − |g(ν₀)−g(ν_h)| → w`,
> so one term is `≥ w/2`. ∎

**Reading.** No continuous (finite-data) instrument reproduces `E` to uniform error `< w/2 = f₁²/2`. To
**certify** `E ≥ r₀` you would need an instrument with error `< r₀/2`, but the instrument's own error is
`≥ f₁²/2 ≥ r₀/2` **at the very point where `E ≥ r₀`**. **The measurement floor equals the signal** —
resolution `≥` signal, circular. This is the exact, quantitative form of ch09 Thm 4 (which is the `w>0`
contrapositive). It closes **model-free** continuous forcing completely.

**The honest caveat (the one exit).** `R3-2`'s supremum is over *all* of `𝓜(K)`, and its proof uses the
bump adversary `ν_h ∉ {phys}`. On the actual **Yang–Mills family** `{ν_{β,L,a}}` the bumps are absent, so
a continuous `g` could in principle match `E` on that thin family **if** one knows enough about it. But
characterizing `{ν_{β,L,a}}` finely enough to build such a `g` **is** the dynamical content — H1/the atom.
That is the corpus's "one exit, two coordinates," now sharpened: **the escape from R3-2 is exactly a
structural assumption on the measure family, and its quantitative force on the atom is the target
itself.** `[O]`.

---

## 3 · The three faces are one obstruction

| face | mechanism | functional that exhibits it | why it is not a proof |
|---|---|---|---|
| **R1 location** | evaluation at `x = m₁² = inf supp ν` | `[C1]` boundary jump, `[C5]` Christoffel, (a) interval projector, (d) resolvent | setting `x` presupposes the gap `m₁²` |
| **R2 limit** | `d,ε → 0` or `n → ∞` | all of `[C1][C3][C5]` | the limit is infinite data = A3; no finite sample reaches it |
| **R3 family** | location-free / structure-free | `[C2]` sup-atom; parametric `[C]` (§5) | targets wrong atom, or assumes low-energy analytic form |

**Theorem R3-1 (edge atom = doubly discontinuous, R1=R2=R3) `[T]`.** `E(ν)=ν({inf supp ν})` is
weak-\* discontinuous at every `ν` with `E(ν)>0`; `E = ev ∘ b` with `b`, `ev` both discontinuous. Any
functional `F` agreeing with `E` on a set containing, for each `ν` with `E(ν)=w>0`, a bump sequence
`ν_h → ν`, is discontinuous by `≥ w/2` (Theorem R3-2) and hence (Lemma 2, density) not a function of any
finite collection of bounded-kernel integrals. Evaluating `F(ν)` requires either the location (**R1**)
or a limit (**R2**) of the exact measure; restricting to avoid the bumps (**R3**) makes `F` continuous
only by **assuming** the family's analytic form — which is the low-energy physics. The three are the same
statement.

---

## 4 · Witness (`research/v3/w7/r3_atom_functional.py`) — quoted, exact numbers

**[C1] Stieltjes boundary, `ν = 0.30 δ₁ + 0.70 δ₃`:** mass-in-window `M(d)` at the atom `x=1` and at an
atomless point `x=2`,

```
   d      M(d)@x=1(atom)   M(d)@x=2(no atom)
  0.5       0.300000          0.000000
 0.01       0.299981          0.000000
0.003       0.299936          0.000000
```

Atom `→ 0.30` exactly; no-atom `→ 0`. **Atomless bump**, mass `0.30`, half-width `h=0.04` at the same
place: `M(0.5)=0.29998, M(0.04)=0.29996, M(0.02)=0.14999, M(0.005)=0.03750` — it drops **only below
`d<h`**. A finite-resolution observer with `d_min > h` sees the **same** `M(d_min)` for atom and bump:
**cannot certify the atom** (`R2`). And the `d=ε` literal formula gives `w/2=0.150`.

**[C2] sup-over-location `A(ν)=sup_x ν({x})`, thin atom `ν_ε=(ε/2)δ_ε+(1−ε/2)δ₂`:** for `ε=10⁻³`,
`A(ν_ε)=0.9995` (the *upper* atom) while the **edge** atom is `ε/2 = 5×10⁻⁴`. Location-free but
targets the **wrong** atom; also `>0` on any finite-volume free field → fails acid (`[NG]`).

**[C3] finite moments blind (`ν_A=0.45δ_{1.0}+0.55δ_{3.3}` vs its atom-smoothed `ν_B`, half-width `h`):**

```
    h      |ΔM0|     |ΔM1|     |ΔM2|     |ΔM3|     |ΔM4|
  0.05    0.0e+00   0.0e+00   8.33e-04  5.66e-03  3.22e-02
 0.001    0.0e+00   0.0e+00   3.33e-07  2.27e-06  1.29e-05
```

`ΔM0 = ΔM1 = 0` **exactly**; `ΔM2 = (A+B)h²/3`, `ΔM3 = (Ax1+Bx2)h²`, … all `→ 0`, yet `f₁²: 0.45 → 0`.
All moments agree while the edge atom jumps by `0.45`. **Stone–Weierstrass:** `sup_ν|E−g| ≥ 0.225 = w/2`
for every weak-\* continuous `g`. `[T]`.

**[C4] thin atom, no finite-data floor:** `ν_ε=(ε/2λ)δ_ε+(1−ε/2λ)δ_{s₂}`, `λ=1`: `χ = w/ε + (1−w)/s₂ =
1.000000` **exactly, for every `ε`** (`ε = 0.5,0.2,0.05,0.01,0.001`), while the edge atom `w=ε/2λ → 0`
and the edge `m = ε → 0`. Any functional of `χ` (the only H1 datum) returns the **same** value for all
`ε` and cannot bound `w` below. **H1 is blind to the atom** — re-derived, matching `H1Sharpening.lean`
(`H1_does_not_imply_atom`).

**[C5] Hankel/Christoffel `λ_n`, `0.30δ_{1/3}+0.70·Unif[0.4,1]`:** `λ_n(atom) = 0.334, 0.303, 0.300,
0.300` (`n=2,4,8,12`) → `w=0.30`; `λ_n(gap) = 0.482, 0.246, 0.0078, 0.0049` → `0`; `λ_n(cont) →`
density. Reads `f₁²` exactly — but only at `x=m₁²` and `n→∞` (**R1+R2**). (At finite `n=4` the *gap*
reads `0.246`, nearly the atom — finite `n` is blind.)

**Solvable acid tests (tool MUST return `0`):**

- **[S1] SU(2) one-plaquette Haar marginal** `ρ₁(z)=π(Cin2π−Cin|z|)` (exact, `ρ₁(0)=7.6581`, `ρ₁'(0)=0`):
  a **smooth continuous density**, no atom. `M(d)@x=0 = 7.59, 3.06, 0.766, 0.153` at `d=0.5,0.2,0.05,0.01`
  → `0` like `2d·ρ₁(0)`. **Correct tool returns `0`** — acid PASS `[V]`. The one-plaquette has *no* edge
  atom (smooth finite-density model).
- **[S2] free massless cut proxy** `ρ(s)=1` on `[0,1]` (atomless): `M(d)@x=0 = 0.2, 0.05, 0.01, 0.003`
  at `d=0.2,0.05,0.01,0.003` → `0` like `d`. **Acid PASS in the continuum `[V]`.** **Caveat `[T]`:** at
  finite `L` the free field is finitely atomic (lowest mode `~(2π/L)²`), `f₁²>0` at finite `L`, `→0`
  only as `L→∞`. **The acid test is an `L→∞` statement** — exactly where the functional is not computable.

---

## 5 · 1D, 2D, and the free field — the solvable-model audit

The solvent models confirm the structure rather than escape it:

- **1D periodic chain `1×L` `[T]` (corpus D6/P6):** finite-range Gibbs measure ⇒ compact transfer ⇒
  **discrete spectrum at finite `L`** ⇒ atoms exist. The boundary/Christoffel reader **returns a positive
  atom** — but that is a **finite-`L` fact** (every finite-volume theory has atoms), not a uniform
  `(L,a)` statement. The gap `λ(L)` is `L`-independent only because 1D is exactly solvable; the reader
  supplies **no** independent certification. `[T]`/`[O]`.
- **2D pure YM `[T-standard]` (Migdal):** exactly solvable, gapped (area law) ⇒ atoms in the continuum.
  The reader returns a positive edge atom, but the positivity is a consequence of **Migdal solvability**,
  a special 2D theorem — not of finite data on `ν`. In 4D this solvability is absent, so the 2D success
  is **not transportable** (the corpus's standing warning: 2D abelian/KT proxies can even have `χ=∞`).
- **Free massless field `[NG]` (the discriminating adversary):** atomless in the continuum (`[S2]`),
  `f₁²=0`; the reader correctly returns `0` — but only as `L→∞`. Any tool "proving" `f₁²>0` here is
  vacuous (acid). This is the model that separates a genuine atom detector from a finite-volume artifact.

**The pattern `[T]`:** on every exactly solvable model the atom reader *works* — because those models are
exactly solvable. On the **one open case (4D `SU(N)`, `β=2.5`)** the reader's evaluation needs the very
gap/isolation it is meant to certify. The solvable models are the **control**, not the proof.

---

## 6 · The one conditional theorem — and its exact circularity `[C]`

The **only** route that produces a *continuous* (finite-data) atom reader is to **restrict the measure
class**. This is worth stating because it pinpoints exactly what "breaking the no-go" costs:

> **Proposition R3-5 (parametric-threshold, `[C]`).** Suppose (A): near `m = inf supp ν`,
> `ν = f·δ_m + c·(s−m)^{θ−1}ds + o((s−m)^{θ−1})` with **`θ ∈ (0,1)` known**, and nothing else in
> `[m, m+I₂)`. Then `(f, m, c)` are determined by the first **three** moments, and `f = E(ν)` is a
> **weak-\* continuous** functional of finitely many moments **on the class (A)**; it is acid-clean
> (`f=0` for the pure power-law cut). So on (A), `f₁² ≥ r₀` is **checkable from finite data**.

**Why this is `[C]`, not a win.** (A) is an **unproven non-perturbative analytic assumption**: it
presupposes the exact low-energy threshold form (a power law with a *possibly present* edge delta and a
clean isolation window `I₂`). The threshold exponent `θ` for 4D `SU(N)` at `β=2.5` is **not known**; the
OPE gives it only perturbatively (the corpus's `FALSIFY-OPE` blindness), and the "nothing else in
`[m,m+I₂)`" clause **is** the isolation/gap structure. So (A)'s content is the low-energy physics, i.e.
the hard part of the target. The proposition shows the **moment of escape** precisely: *the no-go is a
no-go only model-free; the instant you restrict the class enough to make `f` continuously readable, the
restriction is an assumption about the spectrum — circular at the level of the physics, not the
mathematics.* `[C]`, not `[T]`.

---

## 7 · Ledger

| item | statement | label |
|---|---|---|
| edge atom functional | `E(ν)=ν({inf supp ν})=f₁²`, `=ev∘b`, doubly discontinuous | **[T]** |
| discontinuous readers exist | boundary jump `[C1]`, sup-atom `[C2]`, Hankel `[C5]`, blindness `[C3]` | **[T]**+**[V]** |
| atomic acid test | all return `0` on atomless; one-plaquette `[S1]`, free cut `[S2]` | **[V]** |
| **R1 location** | reader needs `x=m₁²=inf supp ν` ⇒ presupposes gap | **[T]** |
| **R2 limit** | reader needs `d,ε→0` or `n→∞` = infinite data = A3 | **[T]** |
| **R3 family** | location-free variant targets wrong atom / assumes analytic form | **[T]**/**[NG]** |
| R1=R2=R3 | one obstruction (edge atom = weak-\* discontinuity) | **[T]** |
| **Stone–Weierstrass bound (Thm R3-2)** | `sup_ν|E−g| ≥ f₁²/2` for every weak-\* continuous `g` | **[T]** |
| ch09 atomweight formula | window `=ε` gives `f₁²/2`; correct = double limit `d↓0` after `ε↓0` | **[T]** (correction) |
| thin atom | `χ=1/λ` for all `ε`, `E(ν_ε)=ε/2λ→0` (H1 blind) | **[T]** (Lean re-check) |
| 1D / 2D | readers positive, but finite-`L`/solvable fact, not uniform | **[T]**/**[O]** |
| free massless | `f₁²=0` in continuum, `>0` at finite `L`; acid is `L→∞` | **[NG]**/**[T]** |
| **computable from exact measure?** | YES (`lim lim`), no gap knowledge needed | **[T]** |
| **computable from finite data?** | NO — SW bound + atomless density | **[NG]** |
| parametric-threshold reader | `[C]`: continuous only by assuming the low-energy form (circular) | **[C]** |
| **H\*** (atom floor + isolation) | irreducible target; readers' uniform value **is** this | **[O]** |
| Clay gap | unbuilt | **[O]** · **0%** |

---

## 8 · Sharpest line

**A discontinuous atom-reading functional exists and is exact — `f₁² = ν({inf supp ν})` is computable from
the full measure by a double limit with no knowledge of the gap, and four analytic machines (Stieltjes
boundary, Christoffel/Hankel, sup-atom, de Branges/spectral-shift) all realize it and all return `0` on
the free field and the one-plaquette. But it is not computable from finite `(L,a)`-data: by the
Stone–Weierstrass resolution bound `sup_ν|E−g| ≥ f₁²/2`, any finite-data instrument has error at least
equal to the signal, so certifying `f₁²≥r₀` is circular by construction (resolution ≥ signal). The three
failure faces — location (presupposes `m₁²`), limit (infinite data = A3), family (assumes the low-energy
analytic form) — are one obstruction, the weak-\* discontinuity of the boundary density at the edge. The
only escape is to restrict to the Yang–Mills family, and the quantitative force of any such restriction on
the atom *is* the target. The exact missing input is unmoved: a uniform-in-`(L,a)` resolution bound below
the isolation `I₂`, i.e. the atom floor with isolation `H*` — which no functional supplies, because every
functional's uniform evaluation is `H*` itself.** **No QED. Clay 0%.**

---

*Refs.* **Atom / weak-\*:** `docs/paper-v2/ch09-weakstar.tex` (Lemmas 1–4, Thm 4, §topology, falsification
table); `research/unconditional/FALSIFY-TOPOLOGICAL.md` (winding = discontinuous ⇒ A3); `ch02-squeeze.tex`
(thin-atom Prop., `eq:squeeze`, `eq:sufficient`, §schur); `lean/YMStatement/H1Sharpening.lean`
(`H1_does_not_imply_atom`, 0 `sorry`). **Tools/classification:** `ch06-tools.tex` (Hardy presence,
Toeplitz/Christoffel upper bound, Nevanlinna, Krein rigidity); `Q6-H1-STATUS.md` (H1 strictly weaker,
carrier hazard); `Q7-TOOL-FORGE.md` (direction theorem, frozen-outlier); `H6-HEAD-SPECTRAL.md`
(typical vs `ess inf`). **Solvable regimes:** `D6-DECISION.md` (1D, 2D, `β=0` Haar); `P6-PERIODIC-TO-H1.md`;
`P1`/`p1_periodic_holonomy.py` (exact SU(2) Haar `ρ₁`). **Standard (textbook, not re-fetched):** Stieltjes–Perron
inversion; Geronimus/Christoffel asymptotics (`liminf λ_n ≥ ν({x})`); Stone–Weierstrass; Birman–Krein /
de Branges spectral-shift and characteristic function; OS reconstruction. **Witness:** `research/v3/w7/r3_atom_functional.py`
(numpy/scipy, one job, no downloads). No `[C]` upgraded to `[T]`; the factor-2 correction is recorded per
charter §2. **No QED; Clay 0%.**
