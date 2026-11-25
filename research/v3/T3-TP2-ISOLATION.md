# T3-TP2-ISOLATION — the order-2 ratio `ρ = m_{0⁺⁺*}/m_{0⁺⁺}` is the **second eigenvalue** of a heat-kernel/Sturm–Liouville transfer operator (exact `ρ²=8/3` in 2D YM, `4` in the rotor, `2` in the TP₂ Mehler kernel); isolation becomes the **gap** only inside a **two-sided window** `1+c ≤ ρ ≤ C`, and the missing half at `β=2.5` is the **lower (TP₂/binding)** bound

> **2026-11-24 · Wave T3 — the order-2 target of the T wave.** Reads `R-WAVE-SYNTHESIS.md`,
> `R2-ISOLATION.md`, `R6-ATOM-ADVERSARIAL.md` (esp. §2.3 double well, §4 R6-a), `R5-EXACT-ANCHORS.md`;
> sibling T-wave witnesses `w7/t6_binding_adversarial.py` (four-cell lattice `A/I/G/B`, double-well
> residue table) and `w7/t1_polymer_rg.py` (Ursell wall), `w7/t2_tube_overlap.py` (RP is order-1);
> `docs/paper-v2/ch02-squeeze.tex`, `ch04-gap.tex`; `research/v3/D6-DECISION.md`,
> `P-WAVE-SYNTHESIS.md`, `P4-PERIODIC-TRANSFER.md`.
> Witness: `research/v3/w7/t3_tp2.py` (`py -3.12 -I`, numpy+scipy, deterministic; 23 checks, all PASS; output §6).
> **Labels (charter §2).** **[T]** theorem/identity · **[C]** conditional · **[NG]** no-go for a named class ·
> **[O]** exact missing step · **[V]** numerical witness. **No QED. Clay 0%.**

---

## 0 · Verdict (up front; no QED)

**The task: compute the order-2 ratio `ρ` exactly in solvable models, name the mechanism, and find a
criterion that could hold at `β=2.5`.** Four statements, each labelled:

1. **The mechanism is the transfer matrix as a heat kernel.** In every solvable confining model the
   0⁺⁺-projected transfer operator is `K_t = e^{−tL}` with `L` a Laplace-type operator, so the
   Källén–Lehmann "masses" are the generator eigenvalues `μ_0 < μ_1 ≤ μ_2 ≤ …`, and
   ```
        ρ² = μ_2 / μ_1 ,     ρ = m_{0⁺⁺*}/m_{0⁺⁺} = √(μ_2/μ_1) .
   ```
   The **order-1** datum (Perron–Frobenius / the first gap) is `μ_1 > 0`; the **order-2** datum
   (isolation) is `μ_2 > μ_1`. Isolation is literally the **second eigenvalue problem**. **[T]**

2. **It is exact and symmetry-fixed in solvable models.** 2D YM (Migdal heat kernel on `SU(N)`):
   `μ_R = (g²/2)C₂(R)`, so `ρ² = C₂(adj)/C₂(fund) = 2N²/(N²−1) ∈ (2, 8/3]` — for `SU(2)` exactly
   `ρ² = 8/3`. The 1D rotor: `μ_n = n²/2I`, `ρ² = 4` in every parity sector. The Mehler (OU) heat
   kernel: `μ_n = n t`, `ρ² = 2`. A particle in a box: `μ_n = n²`, `ρ² = 4`. All are the ratio of two
   eigenvalues of a Laplacian on a compact homogeneous space — no dynamics, no fine tuning. **[T]+[V]**

3. **The structural certificate is total positivity of order 2 (TP₂) — and it is exactly the second
   Karlin rung.** A strictly TP₂ transfer kernel is (unitarily equivalent to) the heat kernel of a
   regular **Sturm–Liouville** operator; Karlin's theorem then gives a simple, positive, **Sturm-ordered**
   spectrum, so `μ_1 < μ_2` is forced. Perron–Frobenius is the order-1 shadow (simple top eigenvalue
   `μ_0 = 0`); TP₂ is what orders `μ_1` against `μ_2`. On the orderless quotient `𝒜/𝒢` there is **no**
   global order, so TP₂ cannot be stated there — the known second-cone no-go. **[T] + [NG]**

4. **The criterion is a two-sided window, and this is the new content.**
   ```
        isolation (fixed margin Δ = μ_2−μ_1 ≥ δ)  +  bounded ratio 1+c ≤ ρ² ≤ C²
        ⇔  a mass gap   μ_1 ≥ δ/(C²−1) > 0 .
   ```
   The **lower** bound `ρ²>1+c` is the TP₂/binding half. The **upper** bound `ρ² ≤ C²` is the
   *level-spacing* half — precisely the `IN-11` heuristic (`excited 0⁺⁺ ≈ 1.5 m_{0⁺⁺}`), now revealed
   as load-bearing rather than cosmetic: the **double well** has `ρ²→∞` (upper bound fails) with the
   isolation margin `Δ = O(1)` fixed and the gap collapsing — *isolation without a gap*. The **free
   field** has no second discrete level at all (`ρ` undefined), so it has neither half. **[T]+[V]**

**Net.** The order-2 ratio is a genuine, exact, model-computable spectral datum; its lower half is a TP₂
statement whose natural home is a *binding-induced* order (the flux tube), and its upper half is a
state-density/`IN-11` statement. **At `β=2.5` the upper half is tractable (no-dissolution, state
density) and the lower half is the `O(u)` binding wall — so the criterion does not reach `β=2.5`.**
It does, however, localize the missing step exactly: a **binding-induced TP₂ second gap** on the coarse
0⁺⁺ kernel. **No QED. Clay 0%.**

---

## 1 · Dictionary: transfer matrix, Perron (order-1) vs TP₂ (order-2)

Let `K_t` be the 0⁺⁺-projected transfer operator over Euclidean time `t`, on the orderless quotient
`𝒜/𝒢`. Write its spectral decomposition `K_t = Σ_n e^{−μ_n t} |n⟩⟨n|` (RP makes `K_t` self-adjoint
and positivity makes the `μ_n` real). The Källén–Lehmann measure of the carrier is `ν = Σ_n w_n δ_{μ_n}`;
the isolation variable is the **generator eigenvalue** `μ_n` (the "mass squared" of the transfer
representation), and the masses are `m_n = √μ_n`.

| rung | object | spectral statement | prototype |
|---|---|---|---|
| **order-1 (Perron)** | the mass gap / first support point | `μ_1 > 0` (simple top eigenvalue `μ_0=0`, positive Perron vector) | PF, Krein–Rutman, Bushell, Bakry–Émery, Cheeger |
| **order-2 (TP₂)** | isolation `I₂` / second support point | `μ_2 > μ_1` (Sturm-ordered simple spectrum) | total positivity of order 2, Karlin flag, Sturm–Liouville |

This is the T3 dictionary: **isolation = the second Perron gap.** Every order-1 tool (PF, RP cone,
Bushell, Bakry–Émery, Cheeger) is blind to `μ_1−μ_2` because `{1,λ₁,λ₁,…}` and `{1,λ₁,λ₂,…}` share every
order-1 datum (R2 §3.3, `PROOF-EDGE-E2`). Seeing `λ₂` requires the **second invariant cone** below the
first — a TP₂ / Sturm structure.

**Two conventions, stated once.** In `d` dimensions the Källén–Lehmann measure lives on `s=m²`; on a
transfer-matrix cylinder the natural spectral variable is the generator eigenvalue `μ=−log λ`, and we
identify `s=μ`, `m=√μ`. In 0+1 (the rotor) the two-point function `⟨O(t)O(0)⟩=Σ w_n e^{−E_n t}` is
already a Laplace transform in the energy gap, so `s=E`. The **ratio `ρ²=μ_2/μ_1`** is convention-free;
its positivity (`>1`) is what the criterion certifies, its exact value is model data.

---

## 2 · Exact computations in solvable models

### 2.1 2D Yang–Mills (Migdal) — the Casimir ratio `[T]+[V]`

The 2D transfer matrix on a spatial circle is the heat kernel on the gauge group, diagonalised by
characters: `K_t = Σ_R e^{−(g²t/2)C₂(R)} Π_R`, with `Π_R` the projector onto the irrep-`R` block. Hence
the generator eigenvalues are
```
   μ_R = (g²/2)·C₂(R) ,     R ranging over the irreps of SU(N).
```
The order-2 ratio is the **Casimir ratio** of the two lightest transfer eigenvalues. In 2D the `R=0`
(trivial) block is one-dimensional — the vacuum only — so the lightest *excitation* is the fundamental
flux loop, the next is the adjoint. Therefore
```
   ρ² = C₂(adj)/C₂(fund) = N / [(N²−1)/(2N)] = 2N²/(N²−1) .
```
`SU(2)`: `C₂(j)=j(j+1)`, `C₂(½)=3/4`, `C₂(1)=2` ⇒ `ρ² = 8/3` exactly. General `N`: `ρ²` runs over
`(2, 8/3]` — **both** window bounds hold with margin (`ρ²>2>1`, `ρ²≤8/3<∞`). `SU(3)`: fund→adj `ρ²=9/4`;
adj→(2,0) `ρ²=10/9` (ratios decrease up the tower — the *first* pair is the binding-sensitive one).

*Carrier honesty* `[T]` (R5, T6 §Part4). `tr F²` in 2D equals the Hamiltonian density `H`, which is
**conserved**, so `⟨0|tr F²|G_R⟩ ∝ C₂(R)⟨0|χ_R⟩ = 0` for `R≠0` by character orthogonality: the `tr F²`
**atom is exactly zero in 2D** (`f₁²=0`). Therefore the `ρ²=8/3` above is the ratio of the **holonomy /
Wilson-loop transfer tower**, *not* the `tr F²` isolation. 2D YM is a clean demonstration of the
**mechanism** (the order-2 ratio is a symmetry-fixed heat-kernel eigenvalue ratio, exactly computable) and
of the **bounded window** — it is *not* a `tr F²`-isolation model. This is the carrier hazard of R6, made
precise; it does not weaken the structural point, but it forbids over-reading the model as 4D evidence.

### 2.2 1D rotor (particle on a ring) — exact tower `[T]+[V]`

`H = −(1/2I)∂_φ²` on `S¹`: `μ_n = n²/(2I)`, `n∈ℤ`. A carrier `O=cos φ` couples `|0⟩` to `n=±1,±2,…`;
`O=cos 2φ` couples only to even `n`. In **both** sectors the ratio of the two lowest *coupled* levels is
```
   ρ² = (2)²/(1)² = 4      (cos φ:  E₂/E₁ = 2/(½)=4)  ;   ρ² = 16/4 = 4   (cos 2φ),
```
i.e. `ρ²=4` exactly, independent of `I`. The free/gapless limit `I→∞` sends every fixed gap `→0`
(`E₁=1/(2I)→0`): a continuum grows from `0`, and the isolation correctly disappears. So the rotor is the
cleanest **positive-and-acid-clean** model: exact bounded window, gap `>0` at finite `I`, no gap in the
free limit.

### 2.3 The TP₂ Mehler heat kernel — the canonical order-2 operator `[T]+[V]`

The Ornstein–Uhlenbeck heat kernel `M_t(x,y) ∝ exp[−((1+q²)(x²+y²)−4qxy)/(2(1−q²))]`, `q=e^{−t}`,
is the textbook **strictly TP₂** kernel: `log M_t` has positive cross-derivative `2q/(1−q²)`, so every
`2×2` minor is `≥0` (log-supermodularity). Its eigenvalues on `L²(dx)` are exactly `λ_n=q^n`, hence
```
   μ_n = −log λ_n = n t ,     ρ² = μ_2/μ_1 = 2   exactly,
```
a **linear tower**. This is the sharpest statement that **TP₂ ⟹ an exact second gap with a fixed ratio**:
TP₂ is the heat kernel of a Sturm–Liouville operator, and Sturm oscillation makes the spectrum simple and
ordered. The anti-correlated Gaussian `exp(−(x²+y²)/2 − bxy)`, `b>0`, has `∂_x∂_y log = −b < 0`,
violates TP₂ (a sampled `2×2` minor is `−0.54`), and is the regime where two levels may cross.

### 2.4 Particle in a box — a *confining* Sturm–Liouville with a bounded window `[T]+[V]`

`−d²/dx²` on `[0,π]` Dirichlet: `μ_n = n²` (`n=1,2,…`), `ρ² = 4`, `μ_1 = 1 > 0`. This is a regular
Sturm–Liouville operator with a *confining* potential; its heat kernel is strictly TP₂; the ratio is
bounded (`C²=4`) and there is a gap. It is the model analogue of "confining ⟹ bounded window + gap" in
the *cleanest* possible setting — and the contrast with the double well (§4) is exactly the difference
between a confining Sturm–Liouville operator and a **double well**, whose two lowest levels tunnel-split.

### 2.5 Model verdict `[T]+[V]`

| model | generator tower | `ρ²` | lower `>1`? | upper `C²<∞`? | gap? |
|---|---|---|---|---|---|
| 2D YM `SU(N)` (holonomy) | `(g²/2)C₂(R)` | `2N²/(N²−1) ∈ (2,8/3]` | ✓ | ✓ | n/a |
| 1D rotor | `n²/2I` | `4` | ✓ | ✓ | ✓ (`I<∞`) |
| Mehler (OU) heat kernel | `n t` | `2` | ✓ | ✓ | ✓ |
| particle in a box | `n²` | `4` | ✓ | ✓ | ✓ |
| double well (Savvidy) | SUSY eigenvalues | `→ ∞` | ✓ | ✗ | ✗ |
| free massless `tr F²` | cut on `[0,∞)` | undefined | ✗ | — | ✗ |

The order-2 ratio is **real, exact, symmetry-fixed, and acid-clean** (absent on the free cut, bounded in
every confining Sturm–Liouville model). The hard content is entirely the 4D continuum.

---

## 3 · The general TP₂ criterion

The structural mechanism is **total positivity of order 2**, and the quantitative criterion is the
Sturm–Liouville gap ratio it entails.

> **Theorem T3-A (TP₂ ⟹ Sturm-ordered spectrum; Karlin/Sturm, `[T-standard]`).** Let `K(x,y)` be a
> continuous, strictly positive kernel on an interval `I` that is strictly TP₂,
> `K(x,y)K(x',y') > K(x,y')K(x',y)` for `x<x'`, `y<y'`. Then `K` is (unitarily equivalent to) the heat
> kernel `e^{−tL}` of a **regular Sturm–Liouville** operator `L=−(p u')'+q u`, `p>0`; its eigenfunctions
> `φ_n` have exactly `n` interior sign changes (Sturm oscillation); the eigenvalues are real, simple and
> ordered, `λ_0>λ_1>λ_2>…>0`, i.e. `μ_0<μ_1<μ_2<…`. In particular `μ_2>μ_1` — **isolation holds**.

This is the sharp content. **Perron–Frobenius (order-1) gives only `λ_0>λ_1`, hence `μ_1>0` — the first
gap. TP₂ (order-2) gives `λ_1>λ_2`, hence `μ_2>μ_1` — the second gap = isolation.** The Mehler (§2.3) and
box (§2.4) are exact realisations; the anti-correlated kernel (§2.3) is a TP₂-violation where the second
gap is unprotected.

Two honest qualifications, both load-bearing:

- **TP₂ is sufficient, not necessary** for `μ_2>μ_1`: a non-TP₂ kernel can still have a simple ordered
  spectrum. TP₂ is the *canonical certificate* and the one that survives coarse graining, not the only
  route.
- **TP₂ gives strictness, not the constant `c`.** The quantitative ratio `μ_2 ≥ (1+c)²μ_1` needs, on top
  of TP₂, a **spectral-gap estimate** of the Sturm–Liouville operator (Cheeger / Lyapunov / a lower bound
  on `q`). In the solvable models the potential is a Laplacian with a known spectrum, so `c` is explicit
  (`ρ² = 8/3, 4, 2`); in 4D it is exactly the open estimate. Stating plain: `TP₂ ⟹ ρ>1` is a theorem;
  `ρ≥1+c` with a uniform `c` is the same `O(u)` cost-vs-rate wall one tier up.

> **Criterion (C-TP2) `[C]`.** Isolation `ρ≥1+c` holds in the 0⁺⁺ channel at coupling `β` if
> **(a) [lower / TP₂]** the coarse 0⁺⁺ transfer kernel is (unitarily equivalent to) the heat kernel of a
> regular Sturm–Liouville operator on an **ordered** one-dimensional *binding coordinate* supplied by the
> flux tube, with first-two-mass ratio `μ_2/μ_1 ≥ (1+c)²`; **and (b) [upper / level spacing]** the ratio is
> bounded, `μ_2/μ_1 ≤ C² < ∞` (the first excited 0⁺⁺ below the two-body cut, `√μ_2 ≤ 2m_g`).
> Then the channel has an isolated atom **and** a gap.

The criterion is **two-sided** — see §4 for why the upper bound cannot be dropped.

---

## 4 · The window theorem — why isolation alone is not the gap `[T]+[V]`

The whole of R2's "SSC equivalence" is one algebraic identity. Write the isolation margin and the ratio:
```
   Δ := Λ★ − m₁² = μ_2 − μ_1 = (ρ² − 1)·μ_1 ,      ρ² = μ_2/μ_1 .
```

> **Theorem T3-B (window equivalence) `[T]`.** Let `ν` be a Källén–Lehmann measure with first atom
> `μ_1>0` (so `m₁=√μ_1>0`) and second support point `μ_2=ρ²μ_1>μ_1` (isolation). Fix `c>0`, `C<∞`, `δ>0`.
> Then
> ```
>   (i)   [1+c ≤ ρ² ≤ C²]  ∧  [Δ ≥ δ]      ⇒      μ_1 ≥ δ/(C² − 1) > 0        (a gap);
>   (ii)  [μ_1 ≥ μ²]  ∧  [ρ² ≥ 1+c]         ⇒      Δ ≥ c·μ² > 0                (isolation);
>   (iii) if ρ² is unbounded, (i) fails: for every δ>0 a measure has Δ=δ and μ_1 = δ/(ρ²−1) → 0  (no gap).
> ```

*Proof.* (i) is `μ_1 = Δ/(ρ²−1) ≥ δ/(C²−1)`. (ii) is `Δ = (ρ²−1)μ_1 ≥ c·μ²` using `ρ²−1 ≥ c` and
`μ_1 ≥ μ²`. (iii) is the same identity read backwards: let `ρ²→∞` at fixed `Δ=δ`; then
`μ_1=δ/(ρ²−1)→0`. ∎

**The double well is the exact witness of (iii)** (witness W4; full residue table in `w7/t6_binding_adversarial.py` §Part2). SUSY
operator `H=−∂²+W`, `W=V'²/4−V''/2`, `V=(g/2)(x²−1)²`:

```
   g      μ_1 = GAP      μ_2        ρ² = μ_2/μ_1     Δ = μ_2−μ_1
  0.5    7.92e−01       3.549           4.5           2.76
  2.0    7.49e−01       6.004           8.0           5.26
  8.0    1.19e−01      13.242         111.2          13.12
 12.0    2.48e−02      20.031         806.3          20.01
 16.0    4.57e−03      27.929        6105.1          27.92
```

As the barrier grows the **gap `μ_1` Kramers-collapses `~e^{−S}`**, the **isolation margin `Δ` stays
`O(1)` and grows**, and the **ratio `ρ²→∞`** — isolation holds in the strongest sense while the gap
vanishes. So:

- **Isolation (lower bound `ρ²>1`, or the absolute margin `Δ>0`) does NOT imply the gap.** The upper
  bound `ρ²≤C²` is *load-bearing*, not cosmetic.
- **`IN-11` (`excited 0⁺⁺ ≈ 1.5 m₁`) is exactly the upper window bound `C≈1.5`.** R2/R6 called it a
  "lattice heuristic"; T3 shows it is the **necessary second conjunct** of the reduced target — without it
  "isolation ⟹ gap" is false (double well). This reframes the R6 verdict "isolation is the gap modulo
  `IN-11`" into an exact statement: *isolation is the gap iff the ratio window is bounded above.*
- **The two acid tests are complementary.** Free field: no second discrete level ⇒ `ρ` **undefined** ⇒
  lower bound fails (W5, acid PASS). Double well: `ρ²→∞` ⇒ upper bound fails ⇒ gap lost (W4). A valid
  criterion must fail or return "undefined" in both — and the window criterion does.

The double well also shows *where* the R6 residue/isolation robustness lives: it is a model with
`f₁²=O(1)`, `Δ=O(1)` **and no gap** (`χ→∞`, H1 fails) — the reason "residue floor" is the tractable half
and "isolation" alone is not Clay.

---

## 5 · No-go on `𝒜/𝒢`, and the additional structure needed

The criterion above cannot be applied to the raw quotient — this is a no-go, and it is what forces the
extra structure.

> **No-go T3-C `[NG]`.** On the orderless quotient `𝒜/𝒢` the coarse transfer kernel carries **no global
> order**, so it cannot be written as a strictly TP₂ / Sturm–Liouville kernel in gauge-quotient
> coordinates; the second invariant cone is **absent** (`INVENT-SECOND-CONE`, `H6-HEAD-SPECTRAL.md` §8).
> Hence criterion (C-TP2)(a) is not available on `𝒜/𝒢`: the lower (order-2) half cannot be *stated*
> there, let alone proved. This is the same "no order on the quotient" wall that blocks FKG / negative
> association / TP₂ in `ch04-gap.tex` §neg-species — isolation and those are one obstruction, tiered.

What supplies the missing order — the task's explicit question, "binding? polymer? reflection
positivity?":

| structure | what it supplies | which half of the window | label |
|---|---|---|---|
| **binding** (R5; `t6_*`) | a one-dimensional **order** (the flux-tube worldline / string modulus) on which the effective operator is Sturm–Liouville | the **lower** (TP₂) half | `[C]` load-bearing |
| **reflection positivity** (`t2_*`) | `K_t` self-adjoint, `K_t=e^{−tL}`, real spectrum, `ν≥0`; the *form* `e^{−tL}` that makes TP₂ meaningful | makes (a) *well-posed* | `[T]` necessary, not sufficient |
| **polymer expansion** (`t1_*`) | reduces the full YM measure to the flux-tube sector with a controlled remainder; controls the quantitative `c` | the **constant** in the lower half | `[C]`, `O(u)` |
| **no-dissolution / state density** | the first excited 0⁺⁺ sits below the two-body cut; `O(1)` window count | the **upper** half | `[C]`/`[V]`, tractable |

**Binding is the discriminant.** RP alone is order-1: log-convexity of `C(t)` holds even when the atom is
`0` (T2 §C), so RP gives the Hilbert space but never the second gap. The polymer expansion converges only
at `β<1/48` (Osterwalder–Seiler; T1's Ursell radius) and is the *transport device*, not the order. Only
**binding** turns the orderless quotient into a problem that *has* an order-2 structure: the flux tube is a
single extended object with a natural longitudinal coordinate, and its worldsheet Hamiltonian is a
Sturm–Liouville (string) operator. Binding reduces the 4D problem to a 1D order-2 problem — which is
precisely why 2D YM and the rotor (whose transfer matrix *is* a heat kernel) are solvable and 4D is not.

---

## 6 · Witness — `research/v3/w7/t3_tp2.py` (`py -3.12 -I`, deterministic; 23 checks, 0 FAIL)

Quoted output (abridged); `numpy 2.x`, `scipy.linalg.eigh_tridiagonal`:

```
W1  2D YM: mu_R=(g^2/2)C2(R).  SU(2) rho^2=2/(3/4)=8/3 EXACT.
    general N: rho^2=2N^2/(N^2-1):  N=2:2.6667  N=3:2.2500  N=100:2.0002  -> 2.
    => rho^2 in (2, 8/3]:  a BOUNDED window (lower>1, upper 8/3<inf).      [T]+[V]
W2  rotor mu_n=n^2/2I:  cos phi -> rho^2=E2/E1=4 ;  cos 2phi -> rho^2=E4/E2=4.
    I->inf: E1=1/(2I)->0  (all gaps -> 0, continuum from 0).                 [T]+[V]
W3  Mehler (OU) heat kernel = strictly TP2 (min 2x2 minor 3.7e-63 >= 0);
    eigenvalues q^n => mu_n = n t  (LINEAR tower) => rho^2 = 2 EXACT.
    anti-correlated Gaussian (b=0.7): min 2x2 minor = -0.54 < 0  (NOT TP2). [T]+[V]
W4  double well (SUSY): gap mu_1 0.79 -> 4.6e-3 (Kramers); Delta=mu_2-mu_1 stays
    O(1) [2.76, 27.92];  rho^2 = mu_2/mu_1 4.5 -> 6105 (upper bound FAILS).
    SSC identity Delta=(rho^2-1)mu_1 to machine precision.                  [T]+[V]
W5  free: massless trF^2 cut from 0 (no atoms, #levels=0 -> rho undefined);
    massive :phi^2: cut from (2m)^2 atomless; elementary phi: 1 level.
    => criterion returns NO isolation.  ACID PASS.                          [T]+[V]
W6  particle in a box -d^2 on [0,pi]: mu_n=n^2, rho^2=4, gap mu_1=1>0
    (confining Sturm-Liouville, BOUNDED window).                            [T]+[V]
W7  window: Delta=1, rho^2<=C^2 => mu_1>=1/(C^2-1): C=1.5 -> 0.80; C=2 -> 0.33.
    C->inf (double well) -> floor->0 (gap not implied).                     [T]
```

All 23 checks PASS. The witness is a check of the inequality structure and the model spectra, not a proof
about Yang–Mills; the theorem content is §3–§4, the numerics are the structural illustrations.

---

## 7 · Reach at `β=2.5` — honest ledger `[O]`

The criterion splits; the two halves have different status at the physical coupling.

- **Upper half (`ρ² ≤ C²`), `IN-11` / level spacing — tractable-ish.** `C` is set by the first excited
  0⁺⁺ relative to the lightest; lattice puts it at `≈1.5–1.7`, below the two-body cut `2m_g`. This is the
  operative-threshold program (`PROOF-ISOLATION §4.4`, `STATE-DENSITY-I2`: `O(1)` window count,
  `L`-uniform) and no-dissolution (`PROOF-NO-DISSOLUTION-KINEMATIC`). It is `[C]`/`[V]`, not proven, but it
  is *the easy half* and it is what fails in the double well.
- **Lower half (`ρ² ≥ 1+c`), TP₂ / binding — the wall.** On `𝒜/𝒢` there is no TP₂ order (§5); the only
  available order is binding-induced, and the quantitative `c` needs the polymer/strong-coupling estimate,
  which converges only at `β<1/48` and cannot cross the Ursell/analyticity wall to `2.5` (T1, R5). This is
  the **same `O(u)` cost-vs-rate residue** that walls the atom, one order-2 tier up (R2 §5, `SUB2`). `[O]`.
- **The criterion is falsifiable at `β=2.5`, though.** R4's minimal GEVP run measures both `m_1` and
  `Δ=m_2−m_1`, hence `ρ` directly. If `ρ` comes out in a bounded window `[1+c, C]` uniformly in `(L,a)`,
  isolation is confirmed *empirically*; a drift of `ρ` toward `1` (level crowding — the dissociation /
  two-copy signature) or toward `∞` (a far level — the double-well signature) would refute the
  isolation-to-gap reading. The `P4` exact 1D transfer matrix can likewise be extended to its **third**
  eigenvalue to test `ρ²=μ_3/μ_2` (the *second* gap) in the tractable chain — a cheap `[V]` check the
  present P4 analysis (top gap `Δ=0.6444`) does not yet report.

**Verdict.** The general criterion is **(C-TP2): a TP₂ (Sturm–Liouville) second gap with a bounded
ratio**, living on a binding-induced 1D order. 2D YM (`ρ²=8/3`) and the rotor (`ρ²=4`) realise it
exactly; the free field and the double well each violate one half and correctly lose the gap. **It does
not reach `β=2.5`** — the lower (TP₂/binding) half is the unbuilt `O(u)` wall; the upper half is in reach
and is exactly `IN-11`. Minimal additional structure needed: **binding** (order), **polymer** (transport
+ the constant `c`), **RP** (the `e^{−tL}` form). Binding is load-bearing.

---

## 8 · Status ledger

| item | statement | label |
|---|---|---|
| order-2 dictionary | isolation `I₂` = second transfer eigenvalue `μ_2>μ_1`; `ρ²=μ_2/μ_1` | **[T]** |
| 2D YM ratio | `μ_R=(g²/2)C₂(R)`; `ρ²=C₂(adj)/C₂(fund)=2N²/(N²−1)∈(2,8/3]`; `SU(2)`:`8/3` | **[T]+[V]** |
| 2D `tr F²` carrier | `tr F²=H` conserved ⇒ atom `f₁²=0`; the `8/3` is the holonomy tower, not `tr F²` | **[T]** |
| rotor ratio | `μ_n=n²/2I`; `ρ²=4` in both parity sectors; free limit gapless | **[T]+[V]** |
| TP₂ ⟹ ordered spectrum | strict TP₂ ⇒ Sturm–Liouville heat kernel ⇒ simple `μ_0<μ_1<μ_2<…` (Karlin) | **[T-standard]** |
| TP₂ is sufficient not necessary | non-TP₂ kernels may still order; TP₂ is the canonical cert | **[T]** |
| `TP₂ ⇒ ρ>1` vs `ρ≥1+c` | strictness is a theorem; the uniform `c` needs a Sturm–Liouville gap estimate | **[T] / [O]** |
| **window theorem** | `Δ=(ρ²−1)μ_1`; isolation + `1+c≤ρ²≤C²` ⇒ gap `μ_1≥δ/(C²−1)` | **[T]** |
| double well = witness | `ρ²→∞` (upper fails), `Δ=O(1)`, gap `→0`; `IN-11` = `C≈1.5` is load-bearing | **[T]+[V]** |
| free field acid test | no second level ⇒ `ρ` undefined ⇒ no isolation | **[T]+[V]** |
| no-go on `𝒜/𝒢` | no global order ⇒ TP₂ absent on the orderless quotient | **[NG]** |
| additional structure | binding (order) + polymer (transport, `c`) + RP (form `e^{−tL}`); binding load-bearing | **[C]** |
| reach at `β=2.5` | upper half (`IN-11`) tractable; lower half (TP₂/binding) = `O(u)` wall | **[O]** |
| Clay gap | unbuilt | **[O] · 0%** |

---

## 9 · Sharpest line

**Isolation `I₂` is the second eigenvalue of the coarse 0⁺⁺ transfer operator, `ρ²=μ_2/μ_1`, and in every
solvable confining model that operator is a heat kernel / Sturm–Liouville operator, so `ρ` is
symmetry-fixed and exact (2D YM `ρ²=8/3`, rotor `4`, Mehler `2`, box `4`, 2D YM general `N`
`2N²/(N²−1)∈(2,8/3]`). The mechanism is total positivity of order 2: TP₂ ⟺ a Sturm–Liouville heat kernel,
whose Sturm oscillation forces `μ_1<μ_2`; Perron–Frobenius is only the order-1 shadow (`μ_1>0`).
But isolation becomes the **gap** only inside a **two-sided window**, `Δ=(ρ²−1)μ_1` with
`1+c≤ρ²≤C²`; the **double well** has `ρ²→∞` and `Δ=O(1)` while the gap collapses — isolation *without* a
gap — so the upper bound, which is exactly the `IN-11` heuristic `excited 0⁺⁺≈1.5m₁`, is **load-bearing**,
not cosmetic. On the orderless quotient `𝒜/𝒢` there is no TP₂ structure at all, so the criterion must run
on a **binding-induced** one-dimensional order (the flux tube), with polymer expansion for the constant
`c` and RP for the `e^{−tL}` form. At `β=2.5` the upper half is in reach and the lower half (binding TP₂)
is the unbuilt `O(u)` wall. **The general criterion is (C-TP2): a binding-induced TP₂ second gap with a
bounded ratio; it is exact in 2D YM and the rotor, fails correctly on the free field and the double well,
and does not reach `β=2.5`. No QED. Clay 0%.**

---

*Refs.* **Order-2 / TP₂:** R2 §3.3 (`PROOF-EDGE-E2`, `INVENT-SECOND-CONE`); `H6-HEAD-SPECTRAL.md` §8 (no
second cone on `𝒜/𝒢`); `docs/paper-v2/ch04-gap.tex` §`sec:negspecies` (TP₂ as a dependence class; no order
on the quotient). **Model spectra:** Migdal (1975) 2D YM (heat kernel on the group); standard
Sturm–Liouville / Karlin total-positivity (Karlin, *Total Positivity*; the OU/Mehler heat kernel); the 1D
rotor and particle-in-a-box (textbook). **Window / SSC:** R2 §5 (`Λ★−m₁²=(ρ²−1)m₁²`, fixed-`δ` `I₂`⟺gap
given bounded `ρ`); R6 §2.3, §4 R6-a, §5 (`IN-11`, double well); `w7/t6_binding_adversarial.py` (four-cell
lattice, double-well table); `PROOF-ISOLATION §4.4` (operative threshold); `STATE-DENSITY-I2.md`;
`PROOF-NO-DISSOLUTION-KINEMATIC.md`. **Binding / transport:** R5 (`PROOF-P3-UNCONDITIONAL`); `w7/t1_polymer_rg.py` (polymer
Ursell wall at `1/48`); `w7/t2_tube_overlap.py` (RP is order-1, log-convex with `f₁²=0`). **First gap / constant mode:** D6
(`χ<∞`⟺saturation, order-1); P4 (exact 1D transfer gap `Δ=0.6444`, the *first* gap). **Witness:**
`research/v3/w7/t3_tp2.py` (23 checks, all PASS, deterministic, `py -3.12 -I`). **No third-party downloads;
no `[C]` upgraded to `[T]`; no QED; Clay 0%.**
