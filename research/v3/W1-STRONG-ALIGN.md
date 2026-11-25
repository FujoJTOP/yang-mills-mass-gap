# W1-STRONG-ALIGN — strong-coupling exact alignment `w₁=s/2`, the next coefficient, and the continuation to `β=2.5`: the next coefficient is **not** available from proven inputs (its `O(u)` part *is* the dressed-amplitude gap), Padé/Borel of the exactly-known normaliser reach `2.5` but are **vacuous**, and `w₁(2.5)>0` survives only conditionally `[C]` on an explicit amplitude headroom, unconditionally only in the large-`N` limit

> **2026-11-24 · Wave W1 (W-WAVE-CHARTER §1).** Task: from U3's exact strong-coupling alignment
> `w₁=s/2`, **compute the next coefficient**, test **Padé / Borel to `β=2.5`**, and state whether `w₁>0`
> **survives**. This note does exactly that, and reports the outcome without embellishment.
> **Reads.** `W-WAVE-CHARTER.md`, `U3-LARGE-N-BINDING.md`, `U4-CONTINUUM-ATOM.md`, `U1-HARNACK-BOUND.md`,
> `T2-TUBE-OVERLAP.md`, `unconditional/PROOF-P3-UNCONDITIONAL.md`, `P3-FOURWALL-DERIVATION.md`,
> `P3-DRESSED-MASS.md`, `P3-FULLREP-KP.md`. **Witness.** `research/v3/w7/w1_strong_coupling.py`
> (`py -3.12 -I`, numpy/scipy + `Fraction`, deterministic, exit 0; output `w1_strong_coupling.out`).
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional/model · **[NG]** no-go
> for a named class · **[O]** exact missing step · **[V]** numerical witness here. **No QED. Clay 0%.**

---

## 0 · Verdict (up front; no QED, Clay 0%)

1. **`w₁(0)=s/2` reproduced `[T]` (mod `[KL]`), with the correct `N`-split.** `w₁=a₀/M₀`,
   `a₀=s/(4N²)`, `M₀(0)=Var((1/N)Re\,tr U)`. For the **complex fundamental (`N≥3`)** `M₀(0)=1/(2N²)`
   ⇒ `w₁(0)=s/2`; for the **real fundamental (`SU(2)`)** `M₀(0)=1/N²` ⇒ `w₁(0)=s/4` (an extra factor
   `2`). U3's "`s/2`" is the `N≥3` statement; the `SU(2)` value is `s/4`. Both are `O(N⁰)`.

2. **The next coefficient does NOT follow from proven inputs `[O]`.** `M₀(β)` is an exact
   one-plaquette moment and gives the **exactly computable** part of the first correction; the
   numerator `f₁²=a₀(1−O(u))` (P3) contributes at the **same order** through the coefficient `c` of
   `O(u)`, and that `c` is **exactly the dressed-amplitude resummation `[R]`** (`P3-DRESSED-MASS`),
   which is the one unprinted step. Concretely:
   ```
   w₁(β) = (s/2)[ 1 − ( 1/2 + c/18 ) β + … ]   (N≥3, SU(3) normaliser)
   w₁(β) = (s/4)[ 1 − ( c/2 ) β + … ]          (SU(2), normaliser even)
   ```
   The **`1/2` (SU(3)) / `0` (SU(2))** pieces are **exact** (from the one-plaquette moment); the
   **`c`** pieces are **[O]**. So "the next coefficient" is **half-exact**: its computable half is
   fixed, its other half is *precisely* the long-standing dressed-amplitude `[R]` gap.

3. **Exact normaliser series (computed here) `[T]`/`[V]`.** `SU(2)` by modified Bessel (exact rational
   coefficients); `SU(3)` by exact Weyl quadrature:
   ```
   SU(2):  M₀(β) = 1/4 − β²/32 + 5β⁴/1536 − 7β⁶/23040 + 13β⁸/491520 − …   (only EVEN powers; d_{abc}=0)
   SU(3):  M₀(β) = 1/18 + β/36 + …                                          (ODD term: d_{abc}≠0)
   ```
   The `SU(3)` linear term is the **cubic Casimir moment** `⟨(Re\,tr U)³⟩_Haar=1/4`; it vanishes for
   `SU(2)` because there is no symmetric rank-3 invariant. This is a **genuine `N`-dependent structure**
   of the first correction that U3's leading `s/2` does not see.

4. **Padé / Borel to `β=2.5` work on the normaliser — and are vacuous `[T]`/`[NG]`.** The
   **one-plaquette normaliser** series has finite radius `|β|≈3.83` (nearest zero of `I₁`), so it
   **converges at `2.5`**: raw partial sums still oscillate at the ~5% level, but Padé `[4/4]` nails
   `M₀(2.5)` to `1e-9`. **But** `M₀(β)=Var(O)>0` for **every** `β`, so the normaliser-only ratio
   `w₁^{norm}=a₀/M₀>0` is **positive by construction at every `β`** — it cannot test survival. A
   continuation that is positive for trivial reasons is **vacuous** (the same failure mode as the
   free-field acid test). This is the sharp `[NG]` for the "continue the normaliser" route.

5. **`β=2.5` is outside the strong-coupling cluster radius at `N=2,3` `[T]`/`[C]`.** The full 4D
   `w₁` is a *cluster* (KP) expansion, whose convergence needs `u·K<1` (`u=c_f/N`, `K≤6.1×10²N`,
   `P3-FULLREP-KP`). At `β=2.5` the exact one-plaquette coefficient gives `u(2.5)=c_f(2.5)/N =
   1.02` (`SU(2)`), `3.55` (`SU(3)`) — **`u≫1`**, and `uK≫1`. For `SU(2)`, `β=2.5` is moreover above
   the known bulk transition `β_c≈2.30`. So the strong-coupling series **cannot be used** at `2.5`
   for physical `N`; only in the **large-`N`** limit, where the leading `u=β/(2N²)=1.25/N²→0`.

6. **Does `w₁>0` survive? — Qualified: `[C]` at `N=2,3`; `[T]`-shaped at large `N`.** `w₁(2.5)>0`
   follows **iff** the dressed-amplitude headroom `C` (in `f₁²=a₀(1−O(u))`, `|c|≤C`) obeys
   `C·u(2.5)<1` **and** the series is inside its radius. At `N=2,3`, `u(2.5)≫1` so the bracket is
   **empty** — it is the **exact missing estimate** (`[O]`), not merely a hard problem. In the
   large-`N` limit `u(2.5)→0`, so for `N>N₀≈760` (crude `K`) the continuation reaches `2.5` and
   `w₁→s/2>0`. This is **not** a Clay statement and **not** uniformity in `N` at fixed `β`.

**One line.** *`w₁(0)=s/2` (`s/4` for `SU(2)`) is exact. The next coefficient splits into an exact
one-plaquette normaliser part (`−1/2` for `SU(3)`, `0` for `SU(2)`) and an **unproven** dressed-amplitude
part that is exactly `P3`'s `[R]` step; Padé/Borel of the exactly-known normaliser reach `β=2.5` but are
positive by construction and therefore vacuous; `w₁(2.5)>0` survives only conditionally `[C]` on an
explicit amplitude headroom `C·u(2.5)<1`, unconditionally only in the large-`N` limit — no QED, Clay 0%.*

---

## 1 · Setup: the exact object and why the "next coefficient" is split

Use the charter/U3/P3 conventions (`U3` §1, `P3-FOURWALL` §2, `T2` §1):
```
S = β Σ_p (1 − (1/N) Re tr U_p) ,   weight ∝ ∏_p e^{(β/N) Re tr U_p} ,
O = (1/N) Re tr U_p ,   M₀(β) = Var_β(O) = ⟨O²⟩_c ,   f₁² = a₀(1−O(u)) ,   a₀ = s/(4N²) ,   u := c_f/N ,
w₁(β) = f₁²/M₀(β) ∈ [0,1] .                                                              (1.1)
```
Here `c_f(β)` is the fundamental character coefficient, `e^{(β/N)Re tr U}=Σ_r c_r χ_r`, `c_f=β/(2N)+
O(β³)` for `N≥3` (doubles to `β/2` for `SU(2)`), and `s∈ℤ_{≥1}` is the `f⊗f̄→1` singlet recoupling
(`P3-FOURWALL` (2.1)). The strong-coupling residue is `[T]` modulo `[KL]` (`P3-CYLINDER-RESIDUE`); the
dressed amplitude `a₀(1−O(u))` is `[C]` on `[R]`.

**The key structural fact.** `M₀(β)` is an **exact one-plaquette moment** — computable in closed form.
`f₁²` carries the **same** `a₀` at leading order but with the `O(u)` amplitude correction whose
coefficient is the unproven dressed-amplitude `[R]`. Since `u=O(β)`, the amplitude correction enters
the **same order** in `β` as the leading normaliser correction. Therefore the first coefficient of
`w₁(β)−s/2` is **not** determined by proven inputs; it splits into an exact part (normaliser) and an
`[O]` part (amplitude).

**Two normalisations, kept distinct — and an `SU(2)` trap.** At Haar,
```
N≥3 (complex fund.):  ⟨|tr U|²⟩=1 ⇒ Var(Re tr U)=(1/4)⟨(χ_f+χ_f̄)²⟩=1/2 ⇒ M₀(0)=1/(2N²) ⇒ w₁(0)=s/2,
SU(2) (real fund.):   tr U real,  ⟨(tr U)²⟩=1           ⇒ Var(tr U)=1          ⇒ M₀(0)=1/N²   ⇒ w₁(0)=s/4.
```
U3's `M₀=1/(2N²)` is the `N≥3` (real and imaginary parts each `1/2`) statement; `SU(2)` gets no `1/2`
because `Re tr U=tr U` is already real. This `SU(2)` factor-2 is **not** a convention slip — it is a
genuine `N`-dependence of `w₁(0)`, and the first `N`-dependent correction below.

---

## 2 · The exact one-plaquette normaliser series (computed) `[T]`/`[V]`

### 2.1 `SU(2)` — exact modified-Bessel coefficients `[T]`

For `SU(2)`, `a_j(β)=(2(2j+1)/β)I_{2j+1}(β)` and `O=(1/2)tr U=(1/2)χ_{1/2}`, `χ_{1/2}²=χ_0+χ_1`, so
```
M₀(β) = (1/4)[ 1 + a₁/a₀ − (a_{1/2}/a₀)² ] ,     ⟨O⟩=(1/2)a_{1/2}/a₀ .                    (2.1)
```
Expanding in `β` (witness `[1]`, `Fraction` arithmetic, exact):
```
M₀(β) = 1/4 − β²/32 + 5β⁴/1536 − 7β⁶/23040 + 13β⁸/491520 − 121β¹⁰/55050240 + 8411β¹²/47563407360 − …
      only EVEN powers  (d_{abc}=0 for SU(2); U↔U† symmetry)                                   (2.2)
M₀(0)=1/4 .
```
The normaliser ratio `g(β)=M₀(0)/M₀(β)=w₁^{norm}(β)/(s/4)`:
```
g(β) = 1 + β²/8 + β⁴/384 − β⁶/11520 + β⁸/737280 + …   (only even terms)                      (2.3)
```

### 2.2 `SU(3)` / generic `N≥3` — the `d_{abc}` linear term `[T]` (SU(3) exact) / `[V]`

For the complex fundamental, `Re tr U=(χ_f+χ_f̄)/2`, `χ_fχ_f̄=1+χ_adj`, and `χ_f²=χ_{[2,0]}+χ_{[1,1]}`.
The exact Haar moments (witness `[2]`, Weyl quadrature, `1e-12`) are
```
⟨(Re tr U)²⟩ = 1/2 ,   ⟨(Re tr U)³⟩ = 1/4 ,   ⟨(Re tr U)⁴⟩ = 3/4 ,   ⟨(Re tr U)⁵⟩ = 15/16 .
```
The **cubic** moment `⟨(Re tr U)³⟩=1/4` is the symmetric rank-3 (`d_{abc}`) Casimir contraction. It gives a
**linear** term in `M₀` (absent for `SU(2)`):
```
M₀'(0) = (1/N²)⟨(Re tr U)³⟩_Haar = 1/36 ,
M₀(β)  = 1/18 + β/36 + O(β²) ,     M₀(0)=1/18 ,   ⇒   g(β)=M₀(0)/M₀(β)=1 − β/2 + O(β²) .   (2.4)
```
(The witness reproduces `M₀(0)=1/18`, `⟨(Re trU)³⟩=1/4`, `M₀'(0)=1/36` to machine precision. The
`SU(3)` series **sign changes** direction later: exact `M₀(2.5)=0.0524 < M₀(0)=0.0556`, so `M₀` is
non-monotone; the linear term is not the whole story.)

**This is the `N`-dependent structure U3's `s/2` does not see.** The first correction to `w₁` beyond
`s/2` differs qualitatively between the real (`SU(2)`, even) and complex (`N≥3`, odd) fundamental.

---

## 3 · The next coefficient: exact normaliser part `+` `[O]` amplitude part

Multiply `f₁²=a₀(1−cu+…)` by `1/M₀(β)=M₀(0)⁻¹(1−m₁β+…)`, `m₁=M₀'(0)/M₀(0)` (and `u=c_f/N≈β/2N²`):
```
w₁(β) = (a₀/M₀(0)) [ 1 − c u ][ 1 − m₁ β + … ]
      = (s/2) [ 1 − c u − m₁ β + … ]   (N≥3) ,     = (s/4) [ 1 − c u + … ]   (SU(2), m₁=0) .
```
With `u≈β/(2N²)`:
```
N≥3 (SU(3)):  m₁ = (1/36)/(1/18) = 1/2 ,  u ≈ β/18  ⇒  w₁(β) = (s/2)[ 1 − ( 1/2 + c/18 ) β + O(β²) ]
SU(2):        m₁ = 0                      ,  u ≈ β/2   ⇒  w₁(β) = (s/4)[ 1 − ( c/2 ) β + O(β²) ] .
                                                                                          (3.1)
```

> **The next coefficient is `−(s/2)(1/2+c/18)` for `SU(3)`, `−(s/4)(c/2)` for `SU(2)`.** The
> **normaliser** pieces `1/2` (`SU(3)`, from the `d_{abc}` cubic moment) and `0` (`SU(2)`) are **exact
> `[T]`** and s-independent. The **`c`** pieces are **exactly the dressed-amplitude `[R]` coefficient**
> of `P3-DRESSED-MASS` — `[O]`. **Proven inputs give only the normaliser half of the next coefficient;
> the other half IS the known `[R]` gap.** There is no further secret: the first correction is
> half-exact and half-missing, and the missing half is already named in the P3 ledger.

**Consistency with §3 of the charter's honesty line.** This is the "precision missing step" the task
asked for — the next coefficient is not computable, *and we can say exactly which quantity is
missing*: the `O(u)` coefficient `c` of the dressed amplitude, equivalently `d f₁²/dβ` at `β=0`.

---

## 4 · Padé / Borel to `β=2.5`: what works, what is vacuous, what is out of radius

Three *different* radius statements, kept separate (a campaign trap — conflating them is exactly how
one manufactures false positivity):

**(A) One-plaquette normaliser series — finite radius `|β|≈3.83`, converges at `2.5` `[T]`/`[V]`.**
`M₀(β)` for `SU(2)` is `(1/4)[1+a₁/a₀−(a_{1/2}/a₀)²]`, a ratio of entire Bessel functions; its Taylor
series has its nearest singularity at the first zero of `I₁(β)` on the imaginary axis, `|β|≈3.8317`.
Hence `β=2.5<3.83` is **inside**. Witness `[4]`: raw partial sums of the even series at `2.5`
oscillate (`0.250 → 0.055 → 0.182 → 0.108 → 0.148 → 0.127` vs exact `0.13412`), but Padé `[2/2]`,
`[3/3]`, `[4/4]` give `0.134138 / 0.134119 / 0.134119` — relative errors `1.4e-4 / 4.6e-7 / 7.8e-10`.
**Padé and Borel genuinely converge on the normaliser at `2.5`.** `[T]`/`[V]`

**(B) But the normaliser-only continuation is VACUOUS `[NG]`.** `M₀(β)=Var_β(O)>0` for **every** `β`
(variance of a non-degenerate variable). Therefore `w₁^{norm}=a₀/M₀>0` **by construction at every `β`**.
A continuation that is positive trivially — independent of any dynamical input — cannot test whether
the **atom** survives; it is the exact analogue of the free-field acid test where `w₁=0` but a naive
positivity tool still "passes." The normaliser is **necessary but cannot certify** `w₁>0`. `[NG]` for
this route.

**(C) The full 4D `w₁` cluster series is OUTSIDE its radius at `β=2.5` for `N=2,3` `[T]`/`[C]`.**
`w₁=f₁²/M₀` is computed by the 4D strong-coupling **cluster (KP)** expansion, whose convergence needs
`u·K<1` (`u=c_f/N`, `K≤80e^{2+β}N≈6.1×10²N`, `P3-FULLREP-KP`). The exact one-plaquette coefficient at
`β=2.5` is `c_f(2.5)=2.042` (`SU(2)`), `10.64` (`SU(3)`), so the expansion parameter is
```
u(2.5)=c_f(2.5)/N = 1.02 (SU(2)) ,  3.55 (SU(3))       ⇒   u ≫ 1 ,  uK ≫ 1 .
```
The `SU(3)` series is visibly divergent at `2.5`: the **linear term alone** `1/18+2.5/36=0.125`
**overshoots** the exact `M₀(2.5)=0.0524` by `140%`. For `SU(2)`, `β=2.5` is above the known bulk
transition `β_c≈2.30` (`[C]`, standard lattice result). **The strong-coupling series cannot be
continued to `2.5` at physical `N`; Padé/Borel built on it are meaningless there.** `[T]`/`[C]`

**Large-`N` is where the continuation IS controlled `[T]` (leading) / `[C]`.** At fixed `β=2.5`, the
**leading** strong-coupling parameter is `u=β/(2N²)=1.25/N²→0` as `N→∞` (deep 't Hooft strong
coupling — `β=2.5` fixed is `λ=2N²/β→∞`). Witness `[6]`: `uK=(1.25/N²)(6.1×10²N)=762.5/N<1` for
`N>763` (crude `K`; the true `K` is smaller, so the true `N₀` is below this). At large `N` the
expansion converges and `w₁(2.5)=s/2+O(1/N²)>0` **given** the amplitude bound (§5). **This is the W2
(large-`N` master field) territory, not a fixed-`N` result.**

---

## 5 · The honest conditional statement and the exact missing estimate

The only obstruction to `w₁(2.5)>0` that remains after the exact normaliser is the dressed-amplitude
coefficient. With `f₁²=a₀(1−cu+…)` and a **proven** bound `|c|≤C`:
```
w₁(2.5) = w₁^{norm}(2.5) · (1 − c·u(2.5)) ,     w₁^{norm}(2.5) = a₀/M₀(2.5) > 0 .
⇒  w₁(2.5) ∈ w₁^{norm}(2.5)·[ 1 − C·u(2.5) , 1 + C·u(2.5) ] ,   lower bound > 0  ⟺  C·u(2.5)<1 .   (5.1)
```
Witness `[5]` (SU(3), `s=1`): `w₁^{norm}(2.5)=0.5299`, `u(2.5)=3.55`; the bracket is positive only if
`C<0.28`. The dressed-amplitude coefficient is `[R]` (unproven); the only explicit constants in hand
(the crude KP entropy `K`, `T2`) are `O(10²)` and give `C·u(2.5)≫1`. **At `N=2,3` the bracket is
empty.** In the large-`N` regime `u(2.5)=1.25/N²≪1`, so any `C=O(1)` yields `C·u<1` and
`w₁(2.5)>0`.

> **The exact missing estimate (precise).** `w₁(2.5)>0` at fixed `N` is equivalent to a
> **provable `O(1)` bound on the dressed-amplitude coefficient `c`** (`P3-DRESSED-MASS` §4 `[R]`),
> **together with** `u(2.5)<1/C` (i.e. inside the cluster radius). Both fail at `N=2,3`, `β=2.5`.
> This is **the same wall** as `T2` §5.1 / `T-WAVE` §5 — W1 does not move it; it **localises** it to one
> coefficient and shows the normaliser is not the obstacle.

**What this is NOT.** It is not a proof of `w₁(2.5)>0` at `N=2,3`; it is not a no-go; and it is not a
statement about the *mass gap* (charter §3: `W` alone proves the atom has a uniform weight, not
isolation or the gap). `γ_{w₁}=0` (`U4`) is independent and unaffected.

---

## 6 · Acid-test battery (charter §2), all five + continuum, with labels

| # | acid test | result here | label |
|---|---|---|---|
| 1 | **free massless field** `w₁=0` | cut-only KL measure, `e^{E₁t}C(t)=2.8e-4→0` (`T=20,40,80`); `a₀=0` (no recoupling); `w₁=0` exactly | **[T]** |
| 2 | **2D YM** `tr F²=H` conserved, `w₁=0` | `⟨χ_{1/2}|χ_{3/2}⟩=1.4e-17` (orthogonality); `f₁²=C₂(R)⟨χ_0|χ_R⟩=0` at every `N`, incl. `N=∞` | **[T]**/**[V]** |
| 3 | **one plaquette** degenerate | no transfer direction: `a₀/M₀=s/2` here is an **artefact**, the "atom" is degenerate (all states one energy); `w₁` is not a spectral fraction — **reported, not invented** | **[T]** |
| 4 | **strong coupling** `w₁=s/2` | `SU(2)`: `s/4` (real fund., `M₀=1/N²`); `N≥3`: `s/2` (complex fund., `M₀=1/(2N²)`) — reproduced exactly | **[T]** (mod `[KL]`) |
| 5 | **large `N`** `w₁=O(N⁰)` | `w₁(0)=s/2` for every `N≥3`, `O(N⁰)`; `u(2.5)→0` as `N→∞` ⇒ continuation controlled at large `N` | **[T]** |
| 6 | **continuum** `γ_{w₁}=0` | `w₁=(Z_O²f₁²)/(Z_O²M₀)` invariant under `O→Z_O O` (`U4`); `dln w₁/dln μ=0` | **[T]** |

All six pass. The `SU(2)` `s/4` vs `s/2` split is the one place the naive reading of "`w₁=s/2`" is
**wrong**; it is corrected here.

---

## 7 · Ledger

| item | statement | label |
|---|---|---|
| `w₁(0)` | `s/2` (`N≥3`, complex fund.), `s/4` (`SU(2)`, real fund.) — `M₀(0)=1/(2N²)` vs `1/N²` | **[T]** (mod `[KL]`) |
| `SU(2)` normaliser | `M₀(β)=1/4−β²/32+5β⁴/1536−…`; `g=M₀(0)/M₀=1+β²/8+β⁴/384−…` (even) | **[T]** |
| `SU(3)` cubic moment | `⟨(Re tr U)³⟩_Haar=1/4` ⇒ `M₀'(0)=1/36`, `g=1−β/2+…` (odd) | **[T]** (SU(3) exact) |
| next coefficient | `−(s/2)(1/2+c/18)` (`SU(3)`), `−(s/4)(c/2)` (`SU(2)`); `1/2`, `0` exact; `c` = `[O]` | **[T]** + **[O]** |
| amplitude `c` | dressed-amplitude `O(u)` coefficient = `P3-DRESSED-MASS` `[R]` gap | **[O]** |
| normaliser Padé/Borel | converge at `2.5` (radius `3.83`; `[4/4]` err `1e-9`) but positive by construction | **[T]**/**[NG]** vacuous |
| full `w₁` at `2.5` | `u(2.5)=1.02` (`SU(2)`), `3.55` (`SU(3)`), `uK≫1`: outside cluster radius; `SU(2)` `β>β_c≈2.30` | **[T]**/`[C]` |
| `w₁(2.5)>0` fixed `N` | `[C]` on `C·u(2.5)<1`; bracket **empty** at `N=2,3` | **[C]**/**[O]** |
| `w₁(2.5)>0` large `N` | `u=1.25/N²→0`; `w₁→s/2>0` for `N>N₀≈760` (crude `K`) | **[T]**-shaped/**[C]** |
| Clay | **unaffected**; W1 localises but does not move the wall (isolation/binding at `2.5`) | — |

---

## 8 · What W1 does and does not establish

- **Does:** reproduce `w₁(0)=s/2` (`s/4` for `SU(2)`) with the exact `N`-split; compute the **exact
  one-plaquette normaliser series** for `SU(2)` and `SU(3)`; exhibit the **cubic-Casimir linear term**
  (`d_{abc}`) that U3's leading result hides; show the next coefficient is **half-exact, half the
  dressed-amplitude `[R]` gap**; demonstrate that Padé/Borel of the normaliser **reach `2.5`** and are
  **vacuous**; and quantify the cluster-radius failure `u(2.5)≫1` at `N=2,3`.
- **Does not:** prove `w₁(2.5)>0` at physical `N`; provide the amplitude bound `C`; or say anything
  about isolation/the gap. The normaliser is **not** the obstacle — the **amplitude coefficient** and
  the **cluster radius** are. The next coefficient's missing half is *exactly* the `[R]` step the
  P3 ledger already names.

*Refs.* `W-WAVE-CHARTER.md`; `U3-LARGE-N-BINDING.md` (§3, §4); `U4-CONTINUUM-ATOM.md` (§2, `γ_{w₁}=0`);
`U1-HARNACK-BOUND.md` (master identity); `T2-TUBE-OVERLAP.md` (§5.1, §6 (F)); internal
`unconditional/PROOF-P3-UNCONDITIONAL.md`, `P3-FOURWALL-DERIVATION.md` (§2), `P3-CYLINDER-RESIDUE.md`,
`P3-DRESSED-MASS.md` (§4 `[R]`), `P3-FULLREP-KP.md` (`K≤6.1×10²N`). External (as verified in P3-AUDIT):
Gupta arXiv:hep-lat/9807028 §14.2.2; Münster PLB 95 (1980)/NPB 190 (1981); Osterwalder–Seiler Ann. Phys.
110 (1978) 440; Seiler LNP 159 (1982); Shen–Zhu–Zhu CMP 400 (2023) 805; Kotecký–Preiss CMP 103 (1986) 491.
Witness `research/v3/w7/w1_strong_coupling.py` (`w1_strong_coupling.out`). **No QED. Clay 0%.**
