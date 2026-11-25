# FOSP-B — `T^{(2)} = ∧²T` is **not** a new positivity structure: it is the **second compound matrix** = TP₂, its spectrum is the pairwise products `{λ_iλ_j}`, and its literal gap is the isolation margin `Δ = μ₂−μ₁` (not the ratio) — **KILL** (reduces to transfer gap + TP₂)

> **2026-11-24 · FOSP-B elimination test.** Reads `T3-TP2-ISOLATION.md`, `W-WAVE-SYNTHESIS.md`,
> `Q6-H1-STATUS.md`. Target: the second exterior power `T^{(2)}=∧²T` as a mechanism to control the
> second spectral layer and force the isolation window `1+c ≤ ρ² ≤ C²`.
> Witness: `research/v3/w7/fosp_b.py` (`py -3.12 -I`, numpy+scipy, deterministic; **15/15 checks PASS**).
> **Labels (charter §2).** **[T]** theorem/identity · **[C]** conditional · **[NG]** no-go for a named class ·
> **[O]** exact missing step · **[V]** numerical witness. **No QED. Clay 0%.**

---

## 0 · Verdict (up front; no QED)

**FOSP-B is KILLED.** `∧²T` is not a new species of positivity; it is the **second compound matrix**, and every
statement it makes about the spectrum is an exact function of `T`'s own spectrum. Four labelled facts:

1. **[T] `∧²T` is the second compound / exactly TP₂.** For a transfer kernel `K(x,y)`, the kernel of `∧²K`
   on antisymmetric functions is the `2×2` minor density `K(x,y)K(x',y')−K(x,y')K(x',y)`. Hence
   `∧²K ≥ 0` **entrywise** ⟺ **all `2×2` minors `≥ 0`** ⟺ **TP₂** (T3 §3). It is the *same* second cone, not a
   third one. (Checks B3a–B3e: Mehler min minor `3.17e−7 ≥ 0` and min `∧²K` entry `3.17e−7`; anti-correlated
   Gaussian min minor `−0.545` and min `∧²K` entry `−0.545` — identical sign.)

2. **[T] No new spectral information.** `spec(∧²T) = {λ_iλ_j}_{i<j}` (B1). The two top eigenvalues are
   `λ₀λ₁=e^{−tμ₁}` and `λ₀λ₂=e^{−tμ₂}`. Its **literal gap** (difference of the two top log-eigenvalues) is
   ```
        Δ_{∧²} = μ₂ − μ₁ = Δ      (the isolation margin, NOT the ratio ρ²),
   ```
   and the ratio recovered from it, `(−log λ₀λ₂)/(−log λ₀λ₁) = μ₂/μ₁`, equals `ρ²(T)` **identically**
   (B2, all `t`, all tested models, `<1e−12`). So `∧²T` carries *no* datum beyond the pair `(μ₁,μ₂)` that `T`
   already carries — it is `T`'s gap and second gap repackaged. **[T]**

3. **[T]+[V] It is blind to the very thing the window needs — the gap.** The window `1+c ≤ ρ² ≤ C²` becomes a
   gap only through `Δ=(ρ²−1)μ₁` **and** the absolute scale `μ₁` (T3-B). `∧²T`'s literal gap is `Δ`, which
   T3-B already proved **insufficient**: the **double well** has `Δ` growing `2.76→27.9` (so "`∧²T` gap"
   looks *healthy*, certifying isolation) while the true gap `μ₁` Kramers-collapses `0.79→4.6e−3` and
   `ρ²=μ₂/μ₁→6111` (B6). `∧²T` would certify isolation in a model with no gap.

4. **[T]+[V] Models: it reproduces T3 exactly and fails exactly where T3 fails.** Rotor `ρ²=4`; Migdal
   SU(2) `ρ²=8/3`; SU(3) fund→adj `ρ²=9/4` (B4) — the *same* numbers as T3 §2, on the *same* (holonomy)
   tower. Free massless circle: `ρ²=4` **independent of `L`** while `μ₁→0` (B5) — the ratio window is
   scale-free and is *satisfied by a gapless free field*; on the `tr F²` atom carrier the spectrum is a
   continuum cut (`f₁²=0`, T3 §2.1, B7) so `ρ` is undefined. Acid pass only because the **gap**, not `∧²T`,
   is what discriminates.

**Reduction.** `∧²T` ⊆ **(transfer gap `μ₁`) + (TP₂ `Δ`)**. The upper window bound `C` (level spacing /
`IN-11` / two-body cut) is still `[O]`; the lower bound's constant `c` still needs the **binding/polymer**
estimate `[O]`; and **on the orderless quotient `𝒜/𝒢` `∧²T` has no entrywise-positive representative**
(T3-C applies verbatim — the no-go is about *order*, which the wedge does not manufacture). **KILL.**

---

## 1 · The identification: `∧²T` is the second compound = TP₂ `[T]`

Let `K_t = e^{−tL} = Σ_n λ_n |n⟩⟨n|`, `λ_n=e^{−tμ_n}`, `μ₀=0<μ₁≤μ₂≤…` (T3 §1). The **second compound**
`∧²K_t` acts on antisymmetric two-vectors with kernel density equal to the `2×2` minor:

```
   (∧²K_t)(x,x'; y,y') = K_t(x,y)K_t(x',y') − K_t(x,y')K_t(x',y) .
```

- **Spectrum.** `spec(∧²K_t) = {λ_iλ_j}_{i<j}` — the second elementary-symmetric (compound) spectrum. This is
  a *function of `spec(T)` alone*: the compound is a functor, `∧²(V⊕W) = ⊕` of the exterior pieces, so no
  new eigenvalue datum is created. **[T]** (B1: six random SPD `T`, max error `<1e−10`; top two compound
  eigenvalues are exactly `λ₀λ₁`, `λ₀λ₂`.)
- **Positivity.** `∧²K_t ≥ 0` entrywise ⟺ every `2×2` minor `≥0` ⟺ **strict TP₂** (T3 §3, Karlin). This is
  literally the statement "the second compound is a positive kernel," i.e. *the* second-cone / TP₂
  certificate — **not** a third, higher rung. **[T]** (B3: Mehler vs anti-correlated Gaussian, identical
  minor/compound signs.)

So the two questions the task poses — "new positivity structure?" and "`∧²T` vs TP₂?" — have the same
answer: **`∧²T` *is* TP₂**, the second compound.

---

## 2 · Why the repackaging cannot close the window `[T]`

The spectral data `∧²T` exposes — `λ₀λ₁=e^{−tμ₁}` (top), `λ₀λ₂=e^{−tμ₂}` (second) — is `(μ₁, μ₂)`.
Two consequences fix the reach:

**(a) The literal `∧²T` gap is `Δ`, the insufficient object.** `−log(top) − (−log(second)) = μ₂−μ₁ = Δ`. T3-B
is the exact statement that `Δ>0` (isolation, which `∧²T` certifies) does **not** imply a gap; the double
well has `Δ=O(1)` rising while `μ₁→0` (B6b). A criterion built on the `∧²T` gap therefore **fails the
double-well acid test** — the same way "isolation alone" does (T3 §4). **[T]+[V]**

**(b) The ratio `∧²T` exposes is identical to `ρ²(T)`.** Reconstructing `ρ² = μ₂/μ₁` from `∧²T`'s top two
eigenvalues returns exactly `ρ²(T)` (B2). So `∧²T` does not *force* `1+c ≤ ρ² ≤ C²` any more than `T` does:
both halves are the same two numbers `(μ₁,μ₂)`, with the same status.

| window half | object | `T` provides | `∧²T` provides | new? |
|---|---|---|---|---|
| lower `ρ²≥1+c` (isolation) | `Δ=μ₂−μ₁` + `c` | `Δ` up to SL/polymer estimate | `Δ` (its literal gap) | **no** (same Δ) |
| upper `ρ²≤C²` (level spacing) | `C` / `IN-11` / two-body cut | nothing analytic | nothing analytic | **no** |
| gap `μ₁` (the scale) | transfer gap | `μ₁` (first gap) | `μ₁` (its top eigenvalue) | **no** |
| constant `c` | binding/polymer `O(u)` | `[O]` | `[O]` | **no** |

Every row is unchanged by passing to `∧²T`. The wedge is a change of basis, not of physics.

---

## 3 · Three elimination tests `[T]+[V]`

### 3.1 Free massless field — the ratio window is *satisfied* by a gapless theory

On a finite circle the free massless transfer tower is `μ_n=(2πn/L)²`. Then `ρ²=μ₂/μ₁=4` **for every `L`**,
while `μ₁=(2π/L)²→0` as `L→∞` (B5). So `1+c ≤ ρ² ≤ C²` **holds** (`4` is in `(2,∞)`) for a theory with no
uniform gap — the ratio is **scale-free** and cannot see `μ₁`. `∧²T` reproduces `ρ²=4` identically. On the
`tr F²` atom carrier the free spectrum is a continuum cut from `0` (`f₁²=0`), so `ρ` is **undefined**
(T3 §2.5). Acid semantics: `∧²T` returns "window OK, gap not certified" — the gap, not the wedge, is the
discriminant. **[T]+[V]** (B5, B7.)

### 3.2 2D holonomy vs `tr F²`

`∧²T` on the 2D Migdal heat kernel gives `ρ²=C₂(adj)/C₂(fund)=2N²/(N²−1)`: SU(2) `8/3`, SU(3) fund→adj
`9/4` (B4) — **exactly** T3 §2.1. But (i) these are the **holonomy/Wilson-loop** tower, and (ii)
`tr F²=H` is conserved in 2D so the `tr F²` atom is `f₁²=0` (T3 §2.1, B7). The `∧²T` construction lives on
the holonomy tower and **cannot reach the `tr F²` isolation channel** — same carrier hazard, unmoved by the
wedge. **[T]+[V]**

### 3.3 SU(3) physical coupling (`β=2.5`)

On the tractable proxies (Migdal exact spectra; the small-volume / strong-coupling transfer tower), `∧²T`
adds no datum: its gap is `Δ=μ_adj−μ_fund`, its ratio is `ρ²=9/4`, both already in `T`'s spectrum. At
`β=2.5` the lower half (`ρ²≥1+c`) is still the **binding TP₂ / polymer `O(u)` wall**, and the upper half
is still **`IN-11` / no-dissolution** — neither is supplied by `∧²T`. **Killed at the physical coupling for
the same reason it is killed in §3.2.** **[V]+[O]**

---

## 4 · Apparent escape hatches, closed

- **"But `∧²T` removes the vacuum — the antisymmetric sector has the second gap directly."** True only as a
  **numerical** device: the Perron mode `|0⟩` is symmetric, so projecting to `∧²` strips it and exposes
  `μ₂−μ₁` with less variance in a GEVP. That is the standard "excited-state / vacuum-removal" trick. It is
  **variance reduction, not a positivity structure**, and it estimates `Δ` (insufficient), not `ρ²` or the
  gap. **[T]**
- **"Positivity of `∧²T` as an operator might be contentful."** No: `∧²` of a positive (self-adjoint,
  positive-spectrum) operator is automatically positive — free of charge, hence **vacuous** (no acid test can
  fail). The only non-trivial positivity is *entrywise* `∧²K ≥ 0`, which is TP₂. **[T]**
- **"A fermionic / Pfaffian RP on `∧²` Fock space is new."** It inherits `RP` (the `e^{−tL}` self-adjoint
  form) and reduces to `RP + TP₂`; it adds no inequality beyond T3's `(C-TP2)`. **[NG]**
- **"Least eigenvalue of `∧²T` bounds `ρ` below/above."** Its spectrum is `{λ_iλ_j}`; any bound on `ρ` is a
  bound on `(μ₁,μ₂)`, i.e. the transfer gap and TP₂ — the named class. **No new functional.** **[NG]**

---

## 5 · Status ledger

| item | statement | label |
|---|---|---|
| identification | `∧²K` kernel density = `2×2` minors; `∧²K≥0` ⟺ TP₂ | **[T]** |
| spectrum | `spec(∧²T)={λ_iλ_j}_{i<j}`; top two `e^{−tμ₁},e^{−tμ₂}` | **[T]** |
| literal gap | `Δ_{∧²}=μ₂−μ₁=Δ` (isolation margin, **not** `ρ²`) | **[T]** |
| no-new-info | `ρ²(∧²T)=ρ²(T)` identically; full spectrum of `∧²T` is a function of `spec(T)` | **[T]** |
| gap blindness | double well: `Δ_{∧²}` rises `2.76→27.9` while `μ₁→0` and `ρ²→6111` | **[T]+[V]** |
| free massless | circle `ρ²=4` ∀`L` while `μ₁→0` (ratio scale-free); `trF²` continuum `ρ` undefined | **[T]+[V]** |
| 2D holonomy vs `trF²` | `∧²T` gives `8/3` (SU(2)), `9/4` (SU(3)) on the holonomy tower; `trF²` atom `f₁²=0` | **[T]+[V]** |
| SU(3) `β=2.5` | no datum beyond `(μ₁,μ₂)`; lower binding/polymer `[O]`, upper `IN-11` `[O]` unchanged | **[V]/[O]** |
| `𝒜/𝒢` | no entrywise-positive `∧²T` representative (T3-C verbatim) | **[NG]** |
| new positivity structure | **does not exist**; `∧²T ⊆` transfer gap + TP₂ | **[NG]** |
| **FOSP-B** | **KILL** (reduces to the transfer gap / TP₂ / H1's building blocks) | **[NG]** |
| Clay gap | unbuilt | **[O] · 0%** |

---

## 6 · Witness — `research/v3/w7/fosp_b.py` (`py -3.12 -I`; 15/15 PASS, deterministic)

```
B1  eig(wedge^2 T) = products of eigs of T; top two = lam0*lam1, lam0*lam2      6 random SPD, err<1e-10
B2  rho^2 from wedge^2 T == rho^2 from T   (rotor, Migdal SU(2)/SU(3), box; t=0.5,1,2.3; err<1e-12)
B3  Mehler: min 2x2 minor 3.17e-7 >=0 ; min entry wedge^2 K 3.17e-7 >=0
    anti-corr Gaussian: min minor -0.545 <0 ; min entry wedge^2 K -0.545 <0   [TP2 <=> wedge^2>=0]
B4  rotor rho^2=4 ; SU(2) 8/3 ; SU(3) fund->adj 9/4      (== T3 sec.2 exactly)
B5  free circle: rho^2=4 for L=4,8,16,32 while mu_1=2.47,0.62,0.15,0.039 -> 0  (ratio scale-free)
B6  double well: mu_1 0.79->4.6e-3 (Kramers); Delta 2.76->27.9 (wedge^2 'gap' healthy); rho^2->6111
B7  trF^2 carrier 2D: f1^2=0 (trF^2=H); wedge^2 lives on holonomy tower, not the atom
TOTAL 15/15 PASS
```

The witness checks the compound identity, the TP₂⟺`∧²≥0` equivalence, the exact model spectra, and the two
acid models (free field, double well). It is a check of the spectral/algebraic structure, not a statement
about Yang–Mills.

---

## 7 · Sharpest line

**`T^{(2)}=∧²T` is the second compound matrix, and that is the whole story: `spec(∧²T)={λ_iλ_j}` is a
function of `spec(T)`, `∧²K≥0` entrywise *is* total positivity of order 2, and its literal gap is the
isolation margin `Δ=μ₂−μ₁` — the exact object T3-B proved insufficient (the double well has `Δ` rising while
the gap dies). `∧²T` reproduces T3's numbers precisely (rotor `4`, SU(2) `8/3`, SU(3) `9/4` on the holonomy
tower, with `trF²` still atomless), certifies isolation but is blind to the gap `μ₁` that the window needs,
gives `ρ²=4` for the gapless free circle, and has no entrywise-positive representative on the orderless
quotient `𝒜/𝒢` (T3-C verbatim). It is not a new positivity species — it is the transfer gap and TP₂ in
disguise, with the antisymmetric projection a numerical vacuum-removal trick. **FOSP-B killed: no new
structure, no reach at `β=2.5`. No QED. Clay 0%.****

---

*Refs.* **T3:** `T3-TP2-ISOLATION.md` (`ρ²=μ₂/μ₁` dictionary §1; 2D Migdal `8/3` §2.1; TP₂⟹Sturm §3;
window theorem T3-B `Δ=(ρ²−1)μ₁` §4; no-go T3-C §5; double-well table §4). **H1:** `Q6-H1-STATUS.md`
(H1 strictly weaker; the gap needs the atom + isolation, not the ceiling). **W:** `W-WAVE-SYNTHESIS.md`
(no theorem at `β=2.5`; `χ_O` sub-linear). **Standard:** exterior powers / compound matrices (Marcus–Minc,
*A Survey of Matrix Theory*); total positivity / Karlin (`T3 §3`); Sturm–Liouville. **Witness:**
`research/v3/w7/fosp_b.py` (15/15 PASS, `py -3.12 -I`). **No third-party downloads; no `[C]` upgraded to
`[T]`; no QED; Clay 0%.**
