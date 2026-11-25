# U3-LARGE-N-BINDING — the atom fraction `w₁=f₁²/M₀` and the isolation ratio `ρ²=μ₂/μ₁` are both **`O(N⁰)`** at large `N`; the `1/N²` single-trace normalisation cancels in the ratio, the multi-trace leakage is `1/N²`-suppressed, binding *does* survive large `N` — but the surviving value is the large-`N` glueball residue / mass ratio `= H1` in the `N`-coordinate, so large-`N` **supports** a uniform binding theorem and is **not** a no-go

> **2026-11-24 · Wave U3 — large-`N` binding.** The T wave left a **nonperturbative binding theorem**: a
> uniform lower bound `w₁≥w₀>0` on the atom fraction of the minimal four-wall tube at `β=2.5`, equivalently
> a leakage bound `Σ_{n≥2}w_n≤1−w₀` (`T-WAVE-SYNTHESIS.md` §5, `T2-TUBE-OVERLAP.md` §8). Large `N` is the
> one limit where the master field / loop equations **close** (`PNT-LARGEN-LOOP`, `Q2-LOOP-BOOTSTRAP`): the
> `1/N` expansion is the solvable-in-principle coordinate. This note computes the **`N`-scaling** of `w₁`
> and of the isolation ratio `ρ²`, tests it on the exact 2D large-`N` theory and on the strong-coupling
> expansion, and decides whether large `N` is evidence for a uniform binding theorem or a no-go.
> **Reads.** `T-WAVE-SYNTHESIS.md`, `T2-TUBE-OVERLAP.md`, `T6-BINDING-ADVERSARIAL.md`,
> `C1-SZZ-NSCALING.md`, `C2-BALABAN-UNIFORM.md`; `unconditional/KAPPA-LARGEN.md`, `PNT-LARGEN-LOOP.md`,
> `CONSTRUCT-MASTERFIELD.md`, `G1-LARGE-N.md`; `v3/Q2-LOOP-BOOTSTRAP.md`, `S5-N-UNIFORM.md`; paper
> `ch05-theorems.tex` (large-`N` / P3), `ch04-gap.tex`. **Witness.** `research/v3/w7/u3_large_n_binding.py`
> (`py -3.12 -I`, numpy only, deterministic, 14/14 PASS; output quoted §8). No third-party downloads.
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional/model · **[NG]** no-go for a
> named class · **[O]** exact missing step · **[V]** numerical witness here. **No QED. Clay 0%.**

---

## 0 · Verdict (up front; no QED, Clay 0%)

1. **`w₁` is `O(N⁰)`, not `N^{-1}` and not `N^{-2}` `[T]+[C]`.** Both the numerator `f₁²` and the
   denominator `M₀=Var(O)` are built from the **same** single-trace operator `O` (e.g. `O=(1/N)tr F²`,
   normalised so `W(C)=O(1)`), so both carry the **same `1/N²`** source normalisation, and the power
   **cancels in the ratio**. Under large-`N` factorization `M₀=Σ_n w_n` is dominated by **single-trace**
   atoms of weight `~g_n/N²` (`O(1)` of them, `g_n=O(1)`), with **multi-trace** atoms at `~1/N⁴` and the
   continuum `~1/N⁴`; hence
   ```
       w₁ = f₁²/M₀ = g₁ / Σ_n g_n + O(1/N²)  =  O(N⁰) ,     f₁² = g₁/N² ,  M₀ = (Σ_n g_n)/N² .
   ```
   The other two candidate scalings are the *individual* objects, not the ratio: `f₁²~N^{-2}` is the
   normalised residue, `N^{-2}` is the **leakage** `Σ_{n≥2}w_n` (multi-trace), `N^{-4}` is `Σ_n w_n` if
   `M₀` were dominated by the bare contact (it is not, after renormalising `O` at a fixed scale). **[T]**
   for the algebraic scaling; **[C]** for the physical applicability (that the large-`N` limit exists and
   `M₀` is exchange-dominated).

2. **`w₁` stays positive as `N→∞` `[C]`.** `w₁∞ = g₁/Σ_n g_n>0` **iff the lightest `0⁺⁺` is a discrete
   single-trace glueball with a nonvanishing residue `g₁>0`**. This is the **generic** large-`N` picture
   (the `0⁺⁺` channel's lightest state is a narrow glueball, width `~1/N→0`; `tr F²` couples to it by
   definition of the channel; no symmetry zeroes `g₁`). So large `N` makes binding the **dominant**
   single-trace contribution, not a fine-tuned one: the atom fraction does **not** degrade. But the
   *value* `w₁∞` is the residue fraction of the lightest glueball — **uncomputed**, because it needs the
   large-`N` spectrum. **[C]** positive · **[O]** value.

3. **The isolation ratio `ρ²` is `O(N⁰)` and bounded `[T]+[C]`.** In the **exactly solvable 2D
   holonomy** carrier `ρ²=C₂(adj)/C₂(fund)=2N²/(N²−1)=2(1+1/(N²−1))→2` **exactly** (SU(2): `8/3`), a
   finite constant `>1` with `|ρ²−2|=2/(N²−1)` — `O(N⁰)`, *bounded away from 1* `[T]`. In 4D `ρ²=m₂²/m₁²`
   is the ratio of the two lightest `0⁺⁺` glueball masses, both `O(Λ)`, `N`-independent, hence `O(N⁰)`
   and bounded `[C]` (conditional on the discrete glueball spectrum). **Large `N` is precisely the
   regime where the window theorem's `ρ²≤C²` holds** (`T3`, `T-WAVE` §1), so isolation `⇒` gap survives.

4. **The large-`N` 2D exact tests pass the acid test `[T]`.** For the **`tr F²` carrier in 2D**,
   `tr F²=tr E²=H` is conserved, `⟨0|tr F²|R⟩=C₂(R)⟨χ₀|χ_R⟩=0` for `R≠0`, so `w₁=0` **exactly at every
   `N`** including `N=∞`: the large-`N` ratio theorem does **not** manufacture an atom — correct, the
   atom is absent because the density is conserved, not because of `N`-counting. For the **holonomy**
   carrier, `ρ²=2N²/(N²−1)` is exact (§3). These are the two exact large-`N` anchors.

5. **Large `N` is NOT solvable in 4D — but the leading-order equation is explicit `[NG]/[O]`.** The
   `N=∞` Makeenko–Migdal equation `∂_μ(δ/δσ_{μν})W=λ∮δ\,W·W` is closed in `W(C)` (no `n`-point tower);
   the `0⁺⁺` residue and mass are the lightest eigenvalue/eigenvector of the **linearised fluctuation
   operator `ℳ` around the master field `W₀`** (`PNT-LARGEN-LOOP` §4). The **missing term** is threefold
   (`Q2` §2.3, §4.4): (i) `W₀` is **unconstructed in 4D** (Singer did 2D only; D=4 surface
   non-uniqueness), (ii) the **discreteness/isolation** of the lightest eigenvalue of `ℳ` is `H1`, and
   (iii) the `O(1/N²)` **connected 4-point `Ĝ`** in (MM) is the multi-trace leakage that shifts `w₁` at
   order `1/N²`. So **large-`N` transports the binding problem faithfully to `N=∞`; it does not reduce it.**

**One line.** *`w₁=f₁²/M₀` and `ρ²=μ₂/μ₁` are both `O(N⁰)`: the `1/N²` single-trace normalisation cancels
in the ratio, so the atom fraction does not degrade with `N` and binding **survives** large `N` (positive
iff the lightest `0⁺⁺` is an isolated single-trace glueball — the generic case); the exact 2D anchors
(`tr F²`: `w₁=0` acid; holonomy: `ρ²→2`) and the strong-coupling `w₁=s/2` confirm it, while the 4D
large-`N` value `w₁∞` is the uncomputed glueball residue/mass `=H1`. **Large `N` is evidence FOR a uniform
binding theorem and NOT a no-go — but not a proof.** No QED. Clay 0%.*

---

## 1 · Set-up: the two ratios, the two normalisations, the acid test

Throughout `O` is a gauge-invariant **single-trace** `0⁺⁺` density — `O=(1/N)tr F²` in the continuum, or
`O=(1/N)Re\,tr U_p` (the plaquette trace) on the lattice — normalised so that the Wilson loop `W(C)=O(1)`
under 't Hooft scaling `λ=g²N` fixed. The Källén–Lehmann measure `ν≥0` on `[0,∞)` is defined by
`C(t)=⟨O(t)O(0)⟩_c=∫e^{−√s\,t}dν(s)`, and (`T2` §1.3, `T6` §1)
```
   f₁² = ν({m₁²}) = |⟨G|O|Ω⟩|² ,   M₀ = Var(O) = ⟨O²⟩_c = ∫dν = Σ_n w_n ,   w_n = |⟨n|O|Ω⟩|² ,
   w₁ = f₁²/M₀ ∈ [0,1]   (the atom fraction; Z_O-invariant, scheme-free — T2 §1.3) ,
   ρ² = μ₂/μ₁ = m₂²/m₁²   (the order-2 isolation ratio — T3; window theorem Δ=(ρ²−1)μ₁, T-WAVE §1) .
```
The task's three candidate scalings are all *formally* possible for `w₁`; §2 shows only `N⁰` is
consistent with the normalisation and with `w₁∈[0,1]`.

**Two normalisations, kept distinct (a campaign trap — `C1-SZZ` §4, `Q2` §1).**
- **SZZ / 't Hooft** `S=Nβ Re\,Σ_p tr U_p`: `K_S=N(1/2−24|β|)`, per-colour constant `N`-free (`C1`).
- **P3 / literal SU(N)** `u\asympβ/N²` (`=β/2N²` generic; `β/N²` SU(2)-special, `P3-FOURWALL` §4): the
  fundamental character coefficient `u` is `1/N²` small at large `N`, so the tube weight `u^{4T}` is
  `N`-small — **but this is the residue, not the ratio** (§4).

**The acid test (charter §1).** A tool that reports a positive `w₁` for the **free massless field** or
the **allowed zero** is vacuous. For large `N` the test is sharpened: the ratio theorem must return
`w₁=0` on the **2D `tr F²` carrier** (conserved density, no atom) — and it does (§3, `A1`), because the
N-counting lives in the *normalisation*, not in the *existence* of the state.

---

## 2 · The core: `w₁` is a ratio of same-normalisation weights, hence `O(N⁰)`

> **Theorem U3-1 (`w₁`, `ρ²` are `O(N⁰)` under large-`N` factorization) `[T]` (algebraic) + `[C]` (applicability).**
> Let `O` be a single-trace `0⁺⁺` operator with the large-`N` normalisation under which `W(C)=O(1)`. (i) The
> connected two-point function factorises over **single-trace** atoms, `M₀=Σ_n w_n`, with single-trace
> weights `w_n≍g_n/N²` (`g_n=O(1)`, `O(1)` of them for `m_n≲2m₁`), **multi-trace** atoms `≍1/N⁴`, and the
> continuum `≍1/N⁴`. (ii) The lightest atom is single-trace, `f₁²=g₁/N²`. (iii) Hence
> ```
>     w₁ = f₁²/M₀ = (g₁/N²) / [(Σ_n g_n)/N² + O(1/N⁴)] = g₁/Σ_n g_n + O(1/N²) = O(N⁰) ,
> ```
> and `ρ²=m₂²/m₁²` is a ratio of two `O(Λ)`, `N`-independent masses, hence `O(N⁰)`. **Both ratios are
> `N⁰` because the `1/N²` source normalisation cancels between numerator and denominator.** ∎

The counting is the standard planar one (it is the same cancellation that makes `K_S∝N`, not `∝N³`, in
`C1` §2, and `κ(N)=κ∞+c/N²` in `KAPPA-LARGEN` §2): a connected correlator of `k` single-trace operators
with the `1/N` normalisation scales as `N^{2−2k}`, so `k=1` (the source `O`) is `O(1)` and `k=2`
(⟨OO⟩_c) is `O(1/N²)`. **The physical `N`-power of `w₁` is set by which weight *dominates* `M₀`, not by
the source.** Three regimes, only the first is large-`N`:

| regime of `M₀` | `f₁²` | `M₀` | `w₁` | meaning |
|---|---|---|---|---|
| **single-trace-dominated** (large `N`) | `g₁/N²` | `(Σg)/N²` | **`g₁/Σg = O(N⁰)`** | the atom is one of `O(1)` single-trace states; **binding survives** |
| multi-trace continuum dominant | `g₁/N²` | `~h/N²`·? | `~N^{-2}` if `M₀` carried an extra `1/N²` from a **fine-tuned light residue** | non-generic (§5 C3) |
| contact-dominated (`M₀~a^{-4}`, bare) | `g₁/N²` | `~a^{-4}` | suppressed by `a` | a UV artifact; removed by renormalising `O` at a fixed scale (`T4` §3) |

**Why not `N^{-1}`.** A `1/N` would require the lightest residue to be one power of `N` *smaller* than
the typical single-trace weight while `M₀` stays single-trace-dominated. There is **no large-`N`
mechanism** for that: all single-trace residues are `≍1/N²` at leading order, and `tr F²` couples to the
`0⁺⁺` glueball with no symmetry zero. So `w₁` is `N⁰`, not `N^{-1}`; the `N^{-1}` is the *amplitude*
`f₁~√(g₁)/N`, not the fraction. **[T]** for the scaling, **[C]** for "no symmetry zero".

---

## 3 · Exact test 1 — 2D large-`N` Yang–Mills (solvable)

2D pure YM is exactly solvable (`Migdal`, `Witten`): `W(C)=e^{−λ\,C₂(\text{fund})\,A(C)/2}` in the
't Hooft normalisation `λ=g²N` (`Q2` §5.1, `PNT-LARGEN-LOOP` §2), and the plaquette variables are
independent. Two exact statements bracket the `tr F²` atom and the isolation ratio.

### 3.1 The `tr F²` carrier: `w₁=0` exactly, all `N` (acid test PASS) `[T]`

In `1+1`D, `tr F²=tr E²=H` is the **energy density**, conserved and diagonal in the character basis.
Under the state–operator map on the cylinder, `⟨0|tr F²|G_R⟩ ∝ C₂(R)\,⟨χ₀|χ_R⟩ = C₂(R)\,δ_{0R}=0` for
`R≠0` (`T2` §4.1, `T6` §4). So the connected `tr F²` two-point function is a `t`-constant **contact**, the
KL measure has **no atom at `m>0`**, and
```
   w₁(tr F², 2D) = f₁²/M₀ = 0   EXACTLY,   for every N including N=∞ .
```
Witness `A1`: SU(2) character orthogonality `∫χ_jχ_k=δ_{jk}` to `6×10⁻¹⁷`; residues `R=0..4` =
`[0.75,0,0,0,0]` (only the vacuum source); `w₁=0`. **This is the decisive acid test for the large-`N`
ratio theorem:** the theorem's `N⁰` scaling is a statement about the *normalisation of the weight*, not
about creating an atom — and correctly returns `0` where the atom is provably absent (conserved density,
2D). The positive recoupling `f⊗f̄→1` is present but yields `w₁=0`, confirming `T2`'s "recoupling is
necessary, not sufficient". **[T]**

### 3.2 The holonomy isolation ratio: `ρ²=2N²/(N²−1)→2` exactly `[T]`

The **holonomy/character** carrier is a *different measure* (`T6` §4, the "2D carrier mismatch"): its
transfer eigenvalues are the quadratic Casimirs `C₂(R)` (the Boltzmann weight `e^{−(g²A/2)C₂(R)}`), so the
isolation ratio between the two lightest nontrivial levels is the **Casimir ratio**
```
   ρ² = C₂(adjoint)/C₂(fundamental) = N / [(N²−1)/(2N)] = 2N²/(N²−1) = 2(1 + 1/(N²−1)) → 2 .
```
Witness `A2`: `ρ²(N)=8/3, 9/4, 32/15, … , 2.000002` at `N=2,3,4,…,1000`; the **exact** deviation
`|ρ²−2|=2/(N²−1)` (to `3×10⁻¹⁶`), bounded in `(2,\,8/3]`. So at large `N`: `ρ²→2`, a **finite constant
`>1`**, `O(N⁰)`, **bounded away from 1**. This is the cleanest exact large-`N` anchor for **isolation**
(not for the `tr F²` atom, which is `0`): **isolation survives large `N` with a bounded ratio** — exactly
the `ρ²≤C²` hypothesis the window theorem needs (§6). **[T]**

---

## 4 · Exact test 2 — strong coupling (P3), `N`-uniformity of `w₁`

At `β<1/48` the four-wall tube residue is explicit (`P3`; `T2` §1.2, §3.1):
```
   f₁² = a₀(1−O(u)) ,   a₀ = (2N)^{-p}\,s > 0 ,   s∈ℤ_{≥1} (singlet recoupling f⊗f̄→1) ,   p∈{0,2} ,
```
and the normaliser `M₀=Var(O)` at `β=0` (Haar one-plaquette) is **exact**:
```
   M₀ = Var((1/N)Re\,tr U) = 1/(2N²)    (uses ∫|tr U|²dU = 1 ⇒ ⟨|tr U|²⟩=1) .
```
Therefore
```
   w₁ = a₀/M₀ = (2N)^{-p}s · 2N² = s·2^{1−p}·N^{2−p} .
```
- `p=0` gives `w₁∝N²`, which **violates `w₁≤1`**: excluded. (Witness `B1`: `w₁=2×10⁸` at `N=10⁴`.)
- `p=2` gives `w₁=s/2 = O(N⁰)`, **`N`-uniform**. (Witness `B2`: `w₁=0.500` for every `N`.)

So the strong-coupling residue floor `w₁≥w₀(β)>0` is **`N`-uniform**: the `(2N)^{-2}` in the amplitude
`a₀` is exactly cancelled by the `1/N²` in `M₀`. This is the **same cancellation** as `C1`'s
`K_S=N(1/2−24|β|)` (Hessian-cost and Ricci-gain both `O(N)`, radius `N`-free) and `KAPPA-LARGEN`'s
`κ(N)=κ∞+O(1/N²)`: the physically meaningful *ratio* is `N`-uniform while the separate *weights* are
`N`-small. The `p=2` matching is **forced** by `w₁≤1`, so it is a *consequence* of the ratio being a
fraction, not an extra input. **[T]** (mod `[KL]` for the spectral identification, as in P3). **The three
candidate scalings are now located exactly:** `w₁~N⁰` (the fraction, `p=2`), `f₁²=a₀~N^{-2}` (the residue),
and `Σ_{n≥2}w_n~N^{-2}` (the leakage) — the last two are the *individual* weights, the first is their
ratio.

---

## 5 · The large-`N` limit in 4D: not solvable, but the leading-order equation and the missing term are explicit

The task asks: if the large-`N` limit is solvable, state the exact result; if not, the leading-order
equation and the missing term. In 4D it is **not** solvable (`PNT-LARGEN-LOOP`, `Q2` §2–§4).

### 5.1 The `N=∞` leading-order equation `[T]`

The large-`N` Makeenko–Migdal equation is **closed in the one functional `W(C)`** (no `n`-point tower —
the precise advantage over the DSE tower, `Q2` §2.2):
```
   ∂_μ (δ/δσ_{μν}(x)) W(C) = λ ∮_C dy_ν δ^{(4)}(x−y) W(C_{xy}) W(C_{yx}) ,      λ = g²N ,        (MM∞)
```
and at finite `N` the `O(1/N²)` **connected 4-point** `Ĝ` is added: `[W·W + (1/N²)Ĝ]` (`Q2` eq. 2.3).
The `0⁺⁺` residue and mass are **not** in `W(C)` itself (that carries the string tension `σ`, the
`1`-form/centre sector) but in the **linearisation of (MM∞) around the vacuum master field `W₀`**
(`PNT-LARGEN-LOOP` §4, `Giles`):
```
   ℳ δW = 0 ,    spec(ℳ) = spec( H_field |_{gauge-invariant} ) ,     m₁² = lightest eigenvalue ,
   f₁²∞ = |⟨G|O|Ω⟩|²  (residue of the lightest pole) ,   M₀∞ = ⟨O O⟩_c(0)  (k=0 limit of ℳ⁻¹) ,
   w₁∞ = f₁²∞ / M₀∞ .
```
So **the exact large-`N` result would be `w₁∞ = g₁/Σ_n g_n` and `ρ²∞ = (m₂/m₁)²`, read off `ℳ`.** In the
rigid-string closure of `CONSTRUCT-MASTERFIELD` these are `m_{0⁺⁺}/√σ=2π√(2κ)` and `ρ²∞` a fixed number —
**provided** the flux-tube rigidity `κ` is known, which is `H1` in string coordinates (`KAPPA-LARGEN` §3).

### 5.2 The missing term(s) `[O]` — three, not one

1. **The master field `W₀` is unconstructed in 4D** (`PNT-LARGEN-LOOP` §3): D=4 spanning surfaces are not
   boundary-determined (a square boundary is spanned by a flat square *and* by tents of every apex
   height), so `W(C)` is a **string functional** (two independent functionals, Coulomb + area), not
   `exp(−σA)`; Singer's construction is 2D only. Without `W₀` there is no `ℳ`.
2. **The isolation/discreteness of `spec(ℳ)` is `H1`** (`Q2` §2.3): even with `W₀`, "is the lightest
   `0⁺⁺` eigenvalue isolated below the next?" is the gap predicate itself, transported to `N=∞`.
3. **The `O(1/N²)` connected 4-point `Ĝ` is the leakage** (`Q2` §4.4): in the curvature coordinate
   `G″=S−D`, `D=Var(force)≥0` is exactly this connected 4-point; bounding it uniformly is `H1`, and using
   large-`N` factorization as a uniform bound is **circular** (`Q2` §4.4; `1/N²=0.25,0.11` at `N=2,3`,
   not small). `Ĝ` is the multi-trace leakage `Σ_{n≥2}w_n` of Theorem U3-1, and it shifts `w₁` only at
   order `1/N²` — but bounding it *uniformly in `(L,a)` at `β=2.5`* is the open input.

**So `w₁∞` is well-defined in the `N=∞` limit (Theorem U3-1: a ratio of `ℳ`-data), but its *value* needs
`ℳ`, which needs `W₀` and the discrete spectrum — i.e. large `N` does not soften the wall, it relocates it.**
This mirrors, exactly, `KAPPA-LARGEN`'s "`κ∞` is the gap in disguise" and `C1`'s "`N`-degradation removed,
the `β`-transport wall remains."

---

## 6 · `ρ²` bounded at large `N` is the window hypothesis: isolation `⇒` gap

The window theorem (`T3`, `T-WAVE` §1): with `Δ=(ρ²−1)μ₁`, isolation together with `1+c≤ρ²≤C²` gives a
uniform floor
```
   μ₁ ≥ δ/(C²−1) ,     δ = the isolation margin (Λ−m₁²) .
```
The double well has `ρ²→∞` (the upper bound `C` is load-bearing) and no floor (`T6` §3). **Large `N`
supplies the missing ingredient structurally:** in 2D the exact holonomy ratio is `ρ²→2` (`C²=2`, §3.2);
in 4D `ρ²→(m₂/m₁)²∞`, a finite `N`-independent constant (the ratio of the two lightest glueball masses),
because both masses are `O(Λ_QCD)` under 't Hooft scaling. So at large `N` **`ρ²=O(N⁰)` and bounded** —
exactly the `ρ²≤C²` side condition under which isolation is a gap. Witness `D1`: `ρ²=2, 8/3` give
floors `Δ/(ρ²−1)=1.0, 0.6>0`; `D2`: the double well `ρ²=4.5→806→…` gives floors `0.79→0.02→0` (the
Kramers collapse). **Reading:** large `N` is the regime in which the *second* half of the binding/isolation
problem (the bounded ratio) is free; the *first* half (the atom exists, `f₁²>0`, and is isolated) is the
large-`N` gap. **[C]** (conditional on the discrete spectrum) for the bounded `ρ²`; **[T]** for the 2D
exact case.

---

## 7 · Does large `N` give evidence for a uniform binding theorem, or a no-go?

**Evidence FOR a uniform binding theorem; NOT a no-go; NOT a proof.** Three positive structural facts, one
honest limit.

**(+) The `N`-scaling is benign in the best possible direction.** `w₁`, `ρ²` are `O(N⁰)`; the leakage is
`N^{-2}` (multi-trace) or `N^{-4}`; the strong-coupling floor is `N`-uniform (`w₁=s/2`); the exact 2D
anchors are finite. This is the *same* benignity `C1` found for SZZ's constant (`λ_LSI=K_S∝N`, per-colour
`N`-free) and `S5` for the block-conditional `(★)` (`δ_W≥cN`, per-colour constant `N`-free). **Large `N`
removes `N`-degradation as an obstruction** — exactly as `C1`/`S5` did for the gap/UV constants. **[T]**

**(+) Binding is *dominant*, not fine-tuned, at large `N`.** The lightest `0⁺⁺` is a narrow single-trace
glueball (width `→0` as `1/N`), the multi-trace cut starts at `2m₁` and is `1/N²`-suppressed, so the atom
fraction is *generic*. The physical intuition (mesons/glueballs are stable narrow states at `N=∞`) and the
scaling agree: **large `N` is where binding is easiest, not where it fails.** **[C]**

**(+) The exact solvable limit is a *positive* calibration.** 2D holonomy `ρ²=2N²/(N²−1)→2` (exact,
bounded, `>1`); 2D `tr F²` `w₁=0` (acid, exactly the vacuous case — the theorem does not over-fire); the
strong-coupling residue `w₁=s/2` (`N`-uniform). The ratio theorem is honest on every solvable case. **[T]**

**(−) The value is the wall.** `w₁∞=g₁/Σg` and `ρ²∞=(m₂/m₁)²` need the large-`N` glueball residue and
mass ratio — the linearised MM spectrum around the **unconstructed** master field, i.e. `H1` in the
`N`-coordinate (`KAPPA-LARGEN` §3, `PNT-LARGEN-LOOP` §4, `Q2` §4.5). Large-`N` factorization is `O(1/N²)`
and **circular** as a uniform bound (`Q2` §4.4: bounding `Ĝ` uniformly *is* `H1`). So the large-`N` route
**does not prove** a uniform binding theorem; it proves the *scaling scaffold* (`N⁰`, positive, leakage
`1/N²`) and relocates the value. **[O]**

**Verdict.** Large `N` is a **scaling theorem, not a proof and not a no-go**: it shows the atom fraction
and the isolation ratio are `O(N⁰)` (binding survives large `N`, dominant and generic), gives exact
solvable anchors, and makes the isolation `ρ²≤C²` side condition **free** at `N=∞`; but its actual value
is the large-`N` glueball residue / mass ratio, i.e. precisely `H1`. **It is genuine evidence *for* a
uniform binding theorem's `N`-structure — the part that says the bound, if it exists, has an `N⁰`
constant — and it is not evidence against (no `N`-counting no-go).** It does not close Clay.

---

## 8 · Witness — `research/v3/w7/u3_large_n_binding.py`

`py -3.12 -I u3_large_n_binding.py` (numpy 2.4.6; exit 0; deterministic; **14/14 PASS**). Quoted output
(abridged):

```
PART A  2D large-N exact
  (A1) SU(2) character orthogonality max|int chi_j chi_k - delta_jk| = 6.18e-17
  (A1) trF^2 residues (R=0..4) = [0.75, 0.0, 0.0, 0.0, 0.0] ; atom fraction w1 = 0.0  (density conserved)
  (A2) rho^2 = C2(adj)/C2(fund) = 2N^2/(N^2-1):
       N=2..1000: 2.666667, 2.250000, 2.133333, ..., 2.000200, 2.000002
  PASS A1a,A1b,A2a,A2b,A2c   (2D trF^2 w1==0 acid; rho^2->2, bounded, |rho^2-2|=2/(N^2-1) exact)
PART B  strong coupling (P3)
  (B1) M0 = Var((1/N)Re tr U) = 1/(2N^2)   (N=2:0.125, N=100:5e-5)
  (B2) w1 = s*2^(1-p)*N^(2-p):  p=0 -> w1~N^2 (excluded, 2e8>1) ;  p=2 -> w1=0.5 (all N)
  PASS B1,B2   (w1<=1 forces p=2; w1 = s/2 = O(N^0), N-uniform)
PART C  large-N Kallen-Lehmann toy (ratio theorem)
   N      w1=f1^2/M0     leakage(MT+ac)/M0     w1(WRONG: f1^2~1/N^4)
   2      0.36036036     3.153e-01            9.009e-02
   1000   0.52631482     1.842e-06            5.263e-07
   1e5    0.52631579     1.842e-10            5.263e-11
  PASS C1 (w1->g1/sum g=0.526316, O(N^0),>0) ; C2 (leakage 1/N^2->0) ;
       C3 (wrong assignment w1~N^-2->0, non-generic) ; C4 (dev O(1/N^2), w1=O(N^0))
PART D  isolation window
  (D1) bounded rho^2: 2.0 -> floor 1.0 ; 8/3 -> 0.6   (>0)
  (D2) double well: g=0.5..12, gap 0.79->0.02, rho^2 4.5->806, floor=gap -> 0.02
  PASS D1,D1,D2   (bounded rho^2 => gap; rho^2->inf => none)
# 14/14 checks PASS
```
Readings: **(A1)** the 2D `tr F²` atom is exactly `0` at all `N` (acid pass — large-`N` counting does
not fake an atom); **(A2)** the holonomy isolation ratio `ρ²→2` exactly, `O(N⁰)`, bounded `>1`; **(B)** the
strong-coupling residue floor is `N`-uniform and `w₁≤1` forces `p=2` ⇒ `w₁=s/2=O(N⁰)`; **(C)** the
large-`N` KL ratio theorem — `w₁→g₁/Σg` (`N⁰`, positive), leakage `1/N²`, the `N^{-2}` outcome only for a
fine-tuned residue; **(D)** bounded `ρ²` (large-`N`) gives a uniform floor, the double well (`ρ²→∞`) does
not. **Witness, not proof; the theorems are §1–§6, the numerics illustrate the scaling and the two exact
2D anchors.**

---

## 9 · Ledger

| item | statement | label |
|---|---|---|
| `w₁` `N`-scaling | `w₁=f₁²/M₀ = g₁/Σ_n g_n + O(1/N²)`, **`O(N⁰)`** (the `1/N²` source normalisation cancels in the ratio) | **[T]** algebraic |
| `f₁²` `N`-scaling | `f₁² ≍ g₁/N²` (normalised single-trace residue); amplitude `f₁~1/N` | **[T]** |
| leakage `Σ_{n≥2}w_n` | multi-trace `≍1/N⁴`, continuum `≍1/N⁴` ⇒ the **`N^{-2}`-type object is the leakage, not `w₁`** | **[T]** |
| `w₁` positivity as `N→∞` | `w₁∞=g₁/Σg>0` iff lightest `0⁺⁺` is an isolated single-trace glueball with `g₁>0` (generic) | **[C]** |
| `ρ²` `N`-scaling | `ρ²=O(N⁰)`, bounded; 4D `=(m₂/m₁)²∞`, `N`-independent | **[C]** |
| 2D `tr F²` atom | `w₁=0` **exactly**, all `N` (`tr F²=H` conserved; `⟨χ₀|χ_R⟩=0`) — acid pass | **[T]** |
| 2D holonomy isolation | `ρ²=C₂(adj)/C₂(fund)=2N²/(N²−1)=2(1+1/(N²−1))→2`; `8/3` at SU(2) | **[T]** |
| strong-coupling `M₀` | `Var((1/N)Re\,tr U)=1/(2N²)` at `β=0` (Haar, `⟨|tr U|²⟩=1`) | **[T]** |
| strong-coupling `w₁` | `a₀=(2N)^{-p}s`, `w₁=a₀/M₀=s·2^{1-p}N^{2-p}`; `w₁≤1⇒p=2⇒w₁=s/2=O(N⁰)`, `N`-uniform | **[T]** (mod `[KL]`) |
| `N=∞` leading equation | (MM∞) closed in `W(C)`; `w₁∞,ρ²∞` from linearisation `ℳ` around master field `W₀` | **[T]** |
| 4D solvability | **not solvable**: `W₀` unconstructed (surface non-uniqueness, Singer 2D only) | **[NG]** |
| missing term(s) | (i) `W₀`/`ℳ`; (ii) discreteness of `spec(ℳ)` `=H1`; (iii) `O(1/N²)` connected 4-pt `Ĝ` (leakage; bounding uniformly is `H1`) | **[O]** |
| window theorem at large `N` | bounded `ρ²=O(N⁰)` is the `ρ²≤C²` hypothesis ⇒ isolation `⇒` gap; 2D exact `ρ²→2`; double well `ρ²→∞` no floor | **[T]** 2D / **[C]** 4D |
| **large-`N` verdict** | evidence **FOR** a uniform binding theorem's `N`-structure (`N⁰` constant, no degradation) · **NOT** a no-go · **NOT** a proof (value `=H1`) | — |
| **Clay** | unbuilt | **[O] · 0%** |

---

## 10 · Sharpest line

**At large `N` the atom fraction and the isolation ratio are `O(N⁰)` — and that is a *structural* fact,
not a dynamical one.** `w₁=f₁²/M₀` is a **ratio of two quantities built from the same single-trace
operator `O`**, so the `1/N²` source normalisation cancels; under large-`N` factorization `M₀` is
single-trace-dominated (`O(1)` atoms of weight `g_n/N²`), the multi-trace leakage is `1/N⁴`, and
`w₁→g₁/Σ_n g_n`, **finite and positive iff the lightest `0⁺⁺` is an isolated single-trace glueball — the
generic large-`N` situation, not a fine-tuned one.** The two exact solvable anchors pin the scaling: 2D
`tr F²` gives `w₁=0` at every `N` (acid: the ratio theorem does not manufacture an atom where the density
is conserved), and the 2D holonomy gives `ρ²=2N²/(N²−1)→2` exactly (bounded, `O(N⁰)`, the `ρ²≤C²` window
hypothesis); the strong-coupling residue `w₁=s/2` is `N`-uniform, the same cancellation that makes SZZ's
constant `∝N` and `κ` `N`-independent. So **binding survives large `N`, the `N`-scaling is benign in the
best direction (`w₁~N⁰`, leakage `~1/N²`), and large `N` supplies the bounded-`ρ²` side condition of the
window theorem for free** — which is genuine evidence **for** a uniform binding theorem and **no**
`N`-counting no-go. But the 4D large-`N` limit is **not** solvable: `w₁∞` is the lightest eigenvalue of the
linearised Makeenko–Migdal operator around the **unconstructed** master field, i.e. the large-`N` glueball
residue and mass ratio `= H1` in the `N`-coordinate, and the `O(1/N²)` connected 4-point that would shift
`w₁` is exactly the leakage whose uniform bound is `H1`. **Large `N` gives the `N`-scaling scaffold of a
uniform binding theorem and its exact 2D calibration; the value remains the gap.** No QED. Clay 0%.

---

*Refs (internal).* `T-WAVE-SYNTHESIS.md` §1,§5; `T2-TUBE-OVERLAP.md` §1.3,§2,§3.1,§4.1,§8;
`T6-BINDING-ADVERSARIAL.md` §1,§3,§4; `T3-TP2-ISOLATION`/`T-WAVE` §1 (window, `ρ²=2N²/(N²−1)`);
`C1-SZZ-NSCALING.md` §2–§4; `C2-BALABAN-UNIFORM.md` §2; `S5-N-UNIFORM.md` §3,§6;
`unconditional/KAPPA-LARGEN.md` §2–§5; `PNT-LARGEN-LOOP.md` §1–§4; `CONSTRUCT-MASTERFIELD.md` §1–§3;
`G1-LARGE-N.md`; `v3/Q2-LOOP-BOOTSTRAP.md` §2.2,§4.4,§4.5; `docs/paper-v2/ch05-theorems.tex`
(P3, `W_geo`, Lemma SC), `ch04-gap.tex` (Savvidy cell, gap predicates).
*External, standard, from the corpus (not re-fetched):* Makeenko–Migdal, Phys. Lett. B88 (1979) 135;
Migdal, Sov. Phys. JETP 42 (1976) 413; Witten, CMP 141 (1991) (2D YM); Singer (2D master field); Giles,
PRD 24 (1981) 2160 (loop spectrum = field spectrum); 't Hooft, NPB 72 (1974) 461 (large-`N` planar
counting). **Witness:** `research/v3/w7/u3_large_n_binding.py` (`py -3.12 -I`, numpy only, deterministic,
14/14 PASS). **No third-party downloads; no `[C]` upgraded to `[T]`; no QED; Clay 0%.**
