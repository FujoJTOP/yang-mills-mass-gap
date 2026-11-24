# C1 — The N-scaling of the Shen–Zhu–Zhu strong-coupling gap / LSI constant

> Attack surface C (LIT-ERA-13 §7 item 3), first angle. Question: *is the N-dependence of the SZZ
> strong-coupling gap/Log-Sobolev constant exponential (large-N route blocked) or polynomial /
> N-independent (publishable lemma)?*
>
> **Verdict (up front).** The N-dependence is **explicit in the paper and polynomial in the
> best possible direction — never exponential.** In the paper's own normalization the Bakry–Émery
> Poincaré/LSI constant is **exactly linear in N**, `K_S = (N+2)/2 − 1 − 8(d−1)N|β|`, positive and
> *increasing* in N, with the strong-coupling radius `|β| < 1/(16(d−1)) = 1/48` **independent of N**
> (the N in the curvature gain and the N in the Hessian cost cancel exactly). The mass-gap
> *correlation-decay rate* that the proof actually extracts is `c_N ≍ K_S/(d_g(a_{e,e}+6(d−1)a_{e,ē}))
> ≍ N^{−3}` — **polynomial degradation**, and the paper says it is "not necessarily optimal."
> So the large-N route is **not blocked by N-degradation of the constant.** The genuine obstruction
> remains the β-transport / continuum-limit wall, which is N-independent.
>
> **Labels:** `[T]` theorem (exact, in the source) · `[T]+[V]` exact + numerically confirmed ·
> `[C]` conditional · `[O]` observation. No `QED` on anything unproven (charter §2).

---

## §1. The source and the extracted theorem

**Citation (VERIFIED).** H. Shen, R. Zhu, X. Zhu, *A stochastic analysis approach to lattice
Yang–Mills at strong coupling*, **Comm. Math. Phys. 400 (2023), no. 2, 805–851**,
DOI 10.1007/s00220-022-04609-1, arXiv:2204.12737 [math.PR]. Confirmed by WebSearch (arXiv abstract,
zbMATH Zbl 1522.81127, INSPIRE 2072947) and against the **arXiv e-print source**
(`Yang_mills_ergodicity_-7.tex`, downloaded 2026-11-24). All formula transcriptions below are from
that source, with line anchors.

**Model (source eq. `measure` / `e:defS`, lines 580–592).** Lattice Yang–Mills on `Λ_L ⊂ ℤ^d`, link
variables `Q_e ∈ G`, `G = SU(N)` or `SO(N)`, Haar measure `σ_N`, action

```
  S(Q) = N β  Re Σ_{p∈P⁺} tr(Q_p),        Q_p = Q_{e1}Q_{e2}Q_{e3}Q_{e4},          (SZZ-1)
```

i.e. the **'t Hooft-scaled** Wilson action (an overall `N` multiplying `Re tr`, with `β` the inverse
't Hooft coupling). The measure is `dμ_{Λ_L,N,β} ∝ e^{S(Q)} ∏_e dσ_N(Q_e)`.

**The ONE hypothesis (source Assumption 1, lines 604–628).**

```
  K_S :=  (N+2)/4 − 1 − 8N|β|(d−1) > 0,      G = SO(N);
          (N+2)/2 − 1 − 8N|β|(d−1) > 0,      G = SU(N).                              (SZZ-2)
```

It is *equivalent* to the strong-coupling condition (source eq. `eq:a`):

```
  |β| <  1/(32(d−1)) − 1/(16N(d−1)),   G = SO(N);
         1/(16(d−1)),                  G = SU(N).        ⇒ d=4:  |β| < 1/48.         (SZZ-3)
```

For **SU(N)** the radius `1/(16(d−1)) = 1/48` in `d=4` is **strictly N-independent** — the `N` in
the positive term `(N+2)/2` and the `N` in the cost `8N|β|(d−1)` are the same power, so `K_S > 0`
does not see `N` to leading order. (SO(N)'s radius does drift with `N`, `→ 1/(32(d−1)) = 1/96` at
`d=4` as `N→∞`; the SU(N) case is the clean one and is the one we quote.)

**Main results (source Theorem 1.3 `th:1.3`, Theorem 1 `th:1`, Theorem 4.2 `th:4.2`).**

| # | Object | Statement (source) |
|---|---|---|
| Thm 1.3 | ergodicity | `W_2^{ρ_{∞,a}}(νP_t, μ_YM) ≤ C(a) e^{−K̃_S t}`; `K̃_S>0` depends on `a,d,β,G` "in particular `N`". |
| Thm 4.2 / eq.(`er`) | Wasserstein gap | `W_2^{ρ_L}(δ_QP_t^L, δ_Q̄P_t^L) ≤ e^{−K_S t} ρ_L(Q,Q̄)`. |
| Thm 1 / eq.(`eqL`) | **LSI** | `μ_YM(F² log F²) ≤ (2/K_S) Σ_{e∈E⁺} μ_YM(|∇_e F|²)`, `μ_YM(F²)=1`. |
| Thm 1 / eq.(`e:PoinYM`) | **Poincaré** | `Var_μ(F) ≤ (1/K_S) Σ_{e∈E⁺} μ_YM(|∇_e F|²)`. |
| Cor. 1 (`co:1`) | large-N Wilson | `Var(W_ℓ/N) ≤ 4n(n−3)/(K_S N)` (SU(N)); `W_ℓ/N → E W_ℓ/N`. |
| Cor. `co:zmm` | **mass gap** | `cov(f,g) ≤ c_1 d_g e^{−c_N d(Λ_f,Λ_g)}(‖f‖_∞‖g‖_∞ + ‖f‖_2‖g‖_2)`, `c_N = c_N(K_S,N,d)`. |

So the paper's own "LSI constant" and "spectral gap" are the **same** object `K_S`, used as a
lower bound (`e^{−K_S t}` / `1/K_S`). It is not hidden, tabulated, or asymptotic — it is **written
out as an equality in Assumption 1.**

---

## §2. The N-counting that produces `K_S` (source lines 1603–1617, 1931–1956)

The whole N-dependence is a two-line computation, and it is the paper's central technical point.

**(a) Cost — the Hessian of the action is only `O(N)`, not `O(N²)`.** Source lines 1604–1611 say,
verbatim in substance: the action "would appear to be of order `N²`" (`tr Q_p ≤ N` times an
overall `N`); *"If the Hessian of S was indeed of order `N²`, or `N^p` for any `p>1`, then … we
would never be able to fix `β` small uniformly in `N` and ensure that `K_S` is strictly positive
when `N` gets large. Fortunately … the Hessian is actually at most of order `N`."* The lemma:

```
  Lemma 4.1 (source eq. e:hess):   |Hess_S(v,v)| ≤ 8(d−1) N |β| |v|²  .
```

This `N¹` (not `N²`) is the **activity/energy count**: `tr D^r(U)` over `r` fundamental factors is
`≤ dim r · unitary`, Hölder arranges it so the `N` from the trace is *not* squared.

**(b) Gain — the Ricci curvature of SU(N) is `O(N)`.** Source eq. `e:Ricvv` (from AGZ (F.6)):

```
  Ric(u,u) = ( α(N+2)/4 − 1 ) |u|²,   α = 1 (SO(N)),  2 (SU(N)).
  ⇒ for SU(N):  Ric(v,v) = ( (N+2)/2 − 1 ) |X|²  =  (N/2) |X|² + O(1).                 (SZZ-4)
```

**(c) Bakry–Émery certification (source eq. `c1`).** The condition `Ric(v,v) − ⟨∇_v∇S,v⟩ ≥ K_S|X|²`
holds with `K_S = [Ric coefficient] − [Hessian coefficient]`:

```
  K_S^{SU(N)} = ( (N+2)/2 − 1 )  −  8(d−1) N |β|  =  N/2 − 8(d−1)N|β| + (1−1)
             = N ( 1/2 − 8(d−1)|β| )  =  N ( 1/2 − 24|β| )      (d = 4).              (SZZ-5)
```

**This is the entire N-dependence. Both the gain and the cost are linear in N; they cancel in the
ratio, so `K_S` itself is linear in N and its positivity radius is N-free.**

---

## §3. The resulting bounds (the answer)

### 3.1 The LSI / Poincaré constant — linear in N (best case) `[T]`

In the task's likely convention, `λ_LSI` is the constant in `Var ≤ (1/λ_LSI) ∫|∇|²` (the
spectral-gap-type *lower bound*; larger = stronger). Comparing with eq. `e:PoinYM`:

```
  λ_LSI(N) = K_S = N ( 1/2 − 24|β| )    for SU(N), d = 4,  |β| < 1/48.                (C1-1)
```

This is a **lower bound of the form `λ_LSI ≥ c·N`** with `c = 1/2 − 24|β| > 0` — strictly stronger
than the polynomial or N-independent regimes the task posed. Equivalently:

- In the `Ent ≤ C_LSI · ∫|∇|²` convention: `C_LSI = 2/K_S = 4/[N(1 − 48|β|)] ∝ 1/N`. Smaller `C`
  = stronger inequality, so **the LSI gets stronger as N grows**.
- **Per color:** `K_S/N = 1/2 − 24|β|` is **exactly N-independent and positive** for `|β|<1/48`.
  This is the clean "N-uniform LSI constant" the campaign wanted.

**No exponential degradation exists.** There is not even polynomial degradation in the LSI/gap
constant — it *improves* linearly.

### 3.2 The mass-gap correlation-decay rate — polynomial `N^{-3}`, not optimal `[T]+[V]`

The paper separates the LSI constant `K_S` from the *mass-gap decay rate* `c_N` of gauge-invariant
observable correlations in lattice spacing (Cor. `co:zmm`). Its remark (source lines 2802–2807)
gives the explicit form:

```
  c_N  ~  K_S / ( d_g ( a_{e,e} + 6(d−1) a_{e,ē} ) ),      "not necessarily optimal",      (C1-2)
```

with (`d_g = dim_ℝ 𝔤 = N²−1`, source Lemma `lem:zmm`, lines 2584–2588, and `γ = 3√2` for SU(N)):

```
  a_{e,ē} = N|β| √d_g ,        a_{e,e} = 2(d−1) N|β| ( √d_g + 6√N ) .                 (C1-3)
```

Substituting `d=4`, `√d_g ≈ N`, `d_g ≈ N²`:

```
  denominator ≈ N² · 6N|β| (4N + 6√N) ≈ 24 N⁴ |β| ,   K_S ≈ N/2 ,
  ⇒  c_N  ≍  (1/2 − 24|β|) / (24 N³ |β|)  =  C(β) / N³ .                              (C1-4)
```

So the **proven** mass-gap decay rate degrades **polynomially as `N^{-3}`** (correlation length
`ξ ≍ N³`). This is the only N-degradation in the paper, and it is an artifact of the crude
commutator estimate `[v_e^i, L]` (`a_{e,·} ∝ N√d_g`, `d_g`), explicitly flagged "not necessarily
optimal." The true lattice mass gap is `≈ 4 ln(1/u)/a`, whose N-scaling is inherited through `u`,
not through this `N^{-3}`.

**Witness (`research/v3/w7/c1_nscaling.py`, run `py -I`).** At `β = 1/96`:

```
   N      K_S(N)     K_S/N      c_N(N)     c_N·N³
   3    7.50e-01   2.50e-01   2.30e-02   6.22e-01
  10    2.50e+00   2.50e-01   6.87e-04   6.87e-01
 100    2.50e+01   2.50e-01   8.70e-07   8.70e-01
  slope d log K_S / d log N = +1.00   (exact: K_S ∝ N)
  slope d log c_N / d log N = −2.90   (→ −3;  c_N·N³ bounded → 0.87)
  c_N(2N)/c_N(N) → 0.126 ≈ 2^{-3}   (polynomial; exponential would → 0)
```

Both halving ratios and the log-log slopes confirm: **K_S grows `∝N`; c_N decays `∝N^{-3}`; neither
is exponential.**

---

## §4. Consistency with P3's `K ≤ CN` and `β₁ = min(1/48, c/N)`

**Same mechanism?** Yes — exactly. `P3-FULLREP-KP.md` §2–§4 obtains, from the full-representation
Kotecký–Preiss hypothesis, an activity bound `|w(P)| ≤ (C_g u)^{|P|}` with `C_g ≤ 2N e^β = O(N)`, and
an animal/entropy factor `λ_a = O(Δ) = O(1)` in N. The KP condition `u·K < 1` then reads
`u · O(N) < 1`. This is the **polymer-expansion image of (SZZ-5)**: the SU(N) group contributes a
factor linear in `N` to the *entropy/gain* (representation-label entropy in P3; Ricci curvature in
SZZ) and a factor linear in `N` to the *activity/cost* (`c_r · dim r` in P3; `|Hess_S|` in SZZ), and
the two are the same power, so **the ratio is N-uniform.** P3's `K ≤ CN` and SZZ's `K_S = O(N)` are
the same counting; the N-cancellation that gives SZZ `β<1/48` is what gives P3 `uK<1` an N-uniform
`u`-radius.

**Does P3 already imply the N-degradation?** **No — it is a coordinate artifact, not a second
mechanism.** `P3-FULLREP-KP.md` §4.3 settles it explicitly: in the *literal* SU(N) normalization the
fundamental coefficient is `u ≍ β/N`, so `uK ≍ (β/N)(CN) = Cβ` and the radius is `β₁ = O(1)` —
**N-independent and consistent with SZZ's `1/48`.** The `c/N` in `β₁ = min(1/48, c/N)` appears only
if one rescales β so that `u ≍ β` ('t Hooft units). Thus the `min(1/48, c/N)` mixes two β
normalizations; neither term is a genuine N-degradation of the *radius*. (There is a real
normalization inconsistency in the campaign's β bookkeeping — SZZ's β multiplies `N·tr`, P3's β
divides `N·tr`, so `β_P3 = N² β_SZZ` — which is exactly what §4.3 of P3-FULLREP-KP is disentangling.
This C1 note does not alter the P3 files; it records that the `c/N` is not an independent wall.)

**The one place P3 *does* carry N-degradation** is the residue *floor* `r₀(β) = ½ c₁ u^p`. Since
`u ≍ β/N` in the literal coordinate, `r₀ ≍ N^{−p}` — **polynomial**, the same structural pattern as
SZZ's `c_N ≍ N^{−3}`: the physically-meaningful *weight* constant degrades polynomially, while the
*radius* (ratio) is N-uniform. Again: polynomial, not exponential.

---

## §5. Acid tests and regime

- **Regime / non-vacuity at N=3.** The statement is strong coupling: `|β| < 1/48 ≈ 0.0208` in `d=4`.
  At `N=3`, `K_S = 3(1/2 − 24|β|) > 0` for `β < 1/48` (e.g. `β=1/96 ⇒ K_S=0.75`), and `c_N > 0`.
  **Positive, non-vacuous constant at N=3**, and for every `N ≥ 2` in the same β window. `[T]`
- **Free-field acid test does not apply — correctly.** The charter's acid test (`f₁²=0`, `G''=0`)
  is for tools that claim to *prove* the gap at weak coupling. SZZ is a strong-coupling
  (high-temperature) theorem and does **not** claim to hold for the free theory at `β→∞`; its
  constant is not asserted there, and indeed `K_S → −8(d−1)N|β|<0` as `β` exits the window. The
  theorem makes no false positivity claim — it is regime-limited, which is honest.
- **The large-N continuum caveat (essential).** SZZ's large-N limit is taken **at fixed `β < 1/48`**
  (deep strong / infinite bare coupling); the paper proves deterministic large-N limits of rescaled
  Wilson loops there. This is the *strong-coupling* 't Hooft limit, **not** the physical continuum
  large-N limit (`λ = 1/β : 48 → 1–3`, `a → 0`). The N-scaling being benign removes N-degradation
  as an obstruction; it does **not** supply the continuum limit. The wall that remains is the
  `β`-transport / `a→0` wall, and it is **N-independent** (same `1/48` for all N).

---

## §6. Verdict and status

| Question | Answer | Status |
|---|---|---|
| Is the N-dependence of the LSI constant stated in SZZ? | **Yes, exactly**, as an equality in Assumption 1: `K_S = (N+2)/2 − 1 − 8(d−1)N|β|`. | `[T]` |
| LSI / Poincaré / spectral-gap constant N-scaling | `λ_LSI = K_S = N(1/2 − 24|β|)`, **linear in N, increasing**; per color `K_S/N = 1/2 − 24|β|` N-independent. | `[T]` |
| Strong-coupling radius N-dependence | `|β| < 1/(16(d−1)) = 1/48`, **N-independent** (SU(N)). | `[T]` |
| Mass-gap decay-rate N-scaling | `c_N ≍ C(β)/N³`, **polynomial**, explicitly "not optimal." | `[T]+[V]` |
| Any **exponential** N-degradation? | **No.** None anywhere in the published constants. | `[T]`+`[V]` |
| Same mechanism as P3's `K ≤ CN`? | **Yes** — activity and entropy both `∝N`, cancel; N-uniform radius. `c/N` in `β₁` is a coordinate artifact. | `[T]`/`[O]` |
| Is the large-N route blocked by N-degradation? | **No.** The constant does not degrade (gap improves, `c_N` polynomial). | `[T]` |
| Does this reach the physical continuum? | **No.** Strong coupling at fixed `β<1/48`; `β`-transport/`a→0` wall remains, N-independent. | `[T]` boundary |

**Publishable lemma (the genuine step).**

> **Lemma C1 (N-uniform strong-coupling LSI; linear-N gap).** For 4D SU(N) lattice Yang–Mills with
> 't Hooft-scaled action `S = Nβ Re Σ_p tr Q_p`, for every `N ≥ 2` and every `|β| < 1/48`
> (the bound independent of N), the Bakry–Émery Poincaré constant is `K_S = N(1/2 − 24|β|) > 0`, so
> the Langevin spectral gap and the LSI constant satisfy `λ_LSI = K_S ≥ (1/2 − 24|β|)·N`; the
> per-color constant `K_S/N = 1/2 − 24|β|` is independent of N. The mass-gap correlation-decay rate
> `c_N` obeys `c_N ≥ C(β)/N³`. In particular the strong-coupling constant does **not** degrade
> exponentially in N, and the large-N strong-coupling limit `N→∞` at fixed `β<1/48` is well-defined
> with a non-vanishing per-color gap.

**Honest boundaries (no QED).** This is a *strong-coupling, fixed-β, fixed-a* statement; it removes
N-degradation as an obstruction but does **not** produce a continuum (`a→0`) or physical-coupling
(`β≈2.5`) mass gap, and the `c_N ≍ N^{-3}` rate is a proof artifact (not physical). Clay is
untouched. `[C]`/`[T]` as labelled; the continuum-mass-gap conclusion is **not** claimed.

---

## §7. Verification ledger

| Item | Check | Status |
|---|---|---|
| Citation identity (authors, title, CMP 400 (2023) 805–851, DOI, arXiv:2204.12737) | WebSearch (arXiv abs, zbMATH Zbl 1522.81127, INSPIRE) **+** read arXiv e-print source | **VERIFIED** |
| Assumption 1 `K_S` formula (SZZ-2) | e-print source lines 604–628 | **VERIFIED (verbatim)** |
| Radius `1/48`, N-independence for SU(N) | source eq. `eq:a` lines 619–628 | **VERIFIED** |
| Hessian bound `|Hess_S| ≤ 8(d−1)N|β||v|²` (Lemma 4.1) | source lines 1613–1616 + the "at most order N" remark lines 1608–1611 | **VERIFIED (verbatim)** |
| Ricci `((N+2)/2 − 1)` for SU(N) (eq. `e:Ricvv`) | source lines 1941–1952 | **VERIFIED** |
| LSI `Ent ≤ (2/K_S)Σμ|∇F|²`, Poincaré `Var ≤ (1/K_S)Σμ|∇F|²` | source Theorem 1, eqs `eqL`,`e:PoinYM`, lines 762–774 | **VERIFIED** |
| Mass-gap rate `c_N ~ K_S/(d_g(a_{e,e}+6(d−1)a_{e,ē}))` | source remark lines 2802–2807; `a_{e,·}` from Lemma `lem:zmm` lines 2584–2588 | **VERIFIED** |
| `K_S ∝ N`, `c_N ≍ N^{-3}`, no exponential | witness `research/v3/w7/c1_nscaling.py` (`py -I`), slopes +1.00 / −2.90, halving ratio → 2^{-3} | **VERIFIED numerically** |

*Primary source: arXiv:2204.12737 e-print `Yang_mills_ergodicity_-7.tex` (isolated extraction in
`research/v3/w7/_szz_src/`). Internal cross-checks: `P3-FULLREP-KP.md` §2–§4, §4.3;
`PROOF-P3-UNCONDITIONAL.md` §3; `PROOF-COARSE-POINCARE.md` §(c) Lemma SC; `LIT-ERA-13.md` entry 13.01
and §7 item 3. No other repository files were edited.*
