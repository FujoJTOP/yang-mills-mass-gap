# W1B-P3-AMPLITUDE — attacking the P3 dressed amplitude `[R]`: the coefficient `c` is **not boundable at `β=2.5`** by any expansion/positivity tool, because it **is the same nonperturbative binding input `W`** — the exact pin `0≤f₁²≤M₀` fixes its window to a factor ≈2 but is **sign-blind**, and the only proven bound on `|c|` (Kotecký–Preiss) is vacuous where `uK≈10⁴`

> **2026-11-24 · Wave W1B (follow-up to `W1-STRONG-ALIGN.md`).** Task: W1 localised the first correction
> of `w₁(β)` to `f₁²=a₀(1−c·u+…)`; the `O(1)` coefficient `c` of the dressed amplitude `D_fin` is exactly
> P3's `[R]` gap. **Can `c` (equivalently the sign/floor of `f₁²`) be bounded at `β=2.5`, or is it the
> same nonperturbative binding input `W`?** Verdict up front. **Reads.** `W1-STRONG-ALIGN.md`,
> `unconditional/PROOF-P3-UNCONDITIONAL.md`, `P3-DRESSED-MASS.md`, `P3-FOURWALL-DERIVATION.md`,
> `T2-TUBE-OVERLAP.md`, `T-WAVE-SYNTHESIS.md`. **Witness.** `research/v3/w7/w1b_p3_amplitude.py`
> (`py -3.12 -I`, numpy/scipy, deterministic, exit 0; output `w1b_p3_amplitude.out`).
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional · **[NG]** no-go for a named
> class · **[O]** exact missing step · **[V]** numerical/one-plaquette witness. **No QED. Clay 0%.**

---

## 0 · Verdict (up front; no QED, Clay 0%)

1. **The object.** `f₁² = A^dressed = a₀·D_fin(β)`, `D_fin := Z_dec(B_T)e^{−f_dec t} = 1 − c·u + O(u²)`,
   `a₀=s/(4N²)`, `u=c_f/N` (`P3-DRESSED-MASS` §3–§4). The `O(1)` coefficient `c` is the **endpoint defect of
   the KP decoration gas** — the one line that P3 proves only as a *bound* `|log D_fin|≤C·u` (P3-DRESSED
   §4.1), valid **iff `u·K<1`**. **It is not a finite-order computable number**: it is the resummed
   decoration self-energy at the tube ends, i.e. intrinsically `[R]`. **[O]**

2. **`[R]` cannot be bounded at `β=2.5`.** The *only* proven bound on `|c|` is the KP absolute-convergence
   bound `|c|≤C`, finite **iff `u·K<1`** with `K(β,N)=80e^{2+β}N` (PROOF-P3 §3). At `β=2.5` the exact
   one-plaquette `c_f` gives
   ```
   u(2.5)·K(2.5) = 1.47×10⁴  (SU(2)) ,   7.66×10⁴  (SU(3))       ⇒   uK ≫ 1 ,
   ```
   while the KP disk closes at `β_KP = 0.0067` (`SU(2)`), `0.0100` (`SU(3)`) — `β=2.5` is `≈250–370×`
   outside. **The convergent constant `C`, and any `O(1)` bound on `c`, lives only on `β≲0.01–0.02`.**
   Witness `[2]`. **[T]`/`[NG]**

3. **The exact window is a factor ≈2 — and sign-blind.** `f₁²` is one atom's weight in the KL measure, so
   `0 ≤ f₁² ≤ M₀ := ∫dμ = Σ_n f_n²` **exactly, at every `β`** (Lemma 2.1). With the exact one-plaquette
   proxy this pins the amplitude at `β=2.5` to
   ```
   f₁²/a₀ ∈ [0, 2.146]  (SU(2)) ,   [0, 1.887]  (SU(3)) .
   ```
   The window **contains both `f₁²=0` and `f₁²≈a₀`** → it cannot decide survival. The two exact anchors
   — strong coupling `f₁²=a₀` at the *midpoint* `w₁=1/2`, and 2D/free `f₁²=0` at the *bottom* — are
   **both consistent with every reflection-positivity / KL / variance / gauge axiom**. Witness `[1]`,`[3]`.
   **[T]**

4. **`[R]` IS `W`. They differ only by the normaliser.** Since `M₀>0`, `sign(w₁)=sign(f₁²)=sign(1−cu+…)`;
   and `w₁=(f₁²/a₀)·(a₀/M₀)` with `a₀/M₀` exact and positive (`SU(3)`, `β=2.5`: `a₀/M₀=0.530`). So
   **`W` ⟺ `f₁²>0`** unconditionally, and a *uniform* amplitude floor `[R]` ⟺ a *uniform* weight floor `W`
   provided the variance is bounded `0<inf M₀≤sup M₀<∞` — and `M₀` is itself the positivity datum that W1
   already showed is **not** the obstruction. Hence `[R]` is the **same nonperturbative binding input** as
   `W` (= `I₂` = `H*` in T2's ledger), not a separate technical gap. Witness `[4]`. **[T]`/`[C]**

5. **Five bound routes fail, none a technicality.** KP radius (§4.1), reflection-positivity/KL positivity
   (§4.2), exact one-plaquette moments (§4.3), absence of an FKG/Griffiths order for `SU(N≥3)` at `β>0`
   (§4.4), and Perron–Frobenius (§4.5) each either dies outside its disk or is **sign-blind** — it returns
   the same value on `f₁²=0` (free field / 2D) and on `f₁²>0` (strong coupling). **[NG]** for the class.

**One line.** *The P3 dressed amplitude `[R]` is the `O(1)` coefficient `c` of a resummed decoration gas,
not a finite-order constant; at `β=2.5` the only proven bound on `|c|` (KP) is vacuous (`uK≈10⁴`), the
exact pin `0≤f₁²≤M₀` fixes the window to a factor ≈2 but is sign-blind (both endpoints realised), and
`[R]` reduces to `W` modulo the exact normaliser — so `[R]` **cannot be bounded at `2.5` by
expansion/positivity tools and is the same nonperturbative binding input**. No QED, Clay 0%.*

---

## 1 · The object: why `c` is `[R]`, not a coefficient one can compute

W1's identity (`W1-STRONG-ALIGN` §3): with `f₁²=a₀(1−cu+…)` and `1/M₀(β)=M₀(0)⁻¹(1−m₁β+…)`,
```
w₁(β) = (s/2)[ 1 − (1/2 + c/18)β + … ]  (SU(3)) ,   (s/4)[ 1 − (c/2)β + … ]  (SU(2)) .        (1.1)
```
The `1/2` (`SU(3)`, from the `d_{abc}` cubic Haar moment `⟨(Re tr U)³⟩=1/4`) and `0` (`SU(2)`) are **exact
one-plaquette** data. The `c` pieces are **not**. By `P3-DRESSED-MASS` §3, `a₀` is the bare four-wall-tube
amplitude (`P3-FOURWALL` §2) and the dressed amplitude is
```
A^dressed = a₀ · D_fin ,  D_fin := Z_dec(B_T) e^{−f_dec t} = 1 + O(u) ,
log D_fin = Σ_{connected Ursell clusters crossing a cut} φ^U(Γ) ,  |log D_fin| ≤ C·u    (P3-DRESSED §4.1).
```
The coefficient `c` in `D_fin=1−cu+…` is the **`O(1)` normalization of the end-cut Ursell sum** — a
transcendental function of the activities, defined only through the resummed decoration-gas pressure. It
is **not** the kind of constant one extracts at finite order: **there is no polynomial in `β` for `c`.** So
the two tasks "compute `c`" and "bound `c`" are the *same* task — control the decoration gas — i.e. `[R]`.
This is a structural, not technical, statement: the amplitude correction *is* the resummation. **[O]**

---

## 2 · The exact pin `0 ≤ f₁² ≤ M₀`, and its sharpness `[T]`

**Lemma 2.1 (two-sided pin). `[T]`** Let `dμ_O≥0` be the KL spectral measure (`P3-CYLINDER-RESIDUE`), so
`M₀ = ∫dμ_O = Σ_n f_n²` and `f₁² = dμ_O({m₁}) = f_{n=1}²`. Since every `f_n²≥0`,
```
0 ≤ f₁² ≤ M₀ = Σ_{n≥1} f_n²    for every β .                                                  (2.1)
```
The left inequality is the atom's positivity, the right is that one term of a sum of nonnegatives is ≤ the
sum. **Both are consequences only of `dμ_O≥0`** — i.e. of reflection positivity and nothing dynamical.

**Both endpoints are attained by exact models `[T]`/`[V]`.** (W1 witness reproduces these.)
```
strong coupling  β=0 :  f₁²=a₀=s/(4N²) ,  M₀=1/(2N²)  (N≥3; =1/N² for SU(2))  ⇒  f₁²/M₀ = s/2 (s/4)  [T]
2D / free        :  f₁²=0 ,  M₀=χ_O>0      (tr F²=H conserved, ⟨χ_{1/2}|χ_{3/2}⟩=1.4×10⁻¹⁷)         [T]/[V]
```
At `β=0` the strong-coupling atom sits at the **midpoint** (`w₁=1/2`) of `[0,M₀]`; the 2D/free absence at
the **bottom** (`f₁²=0`). Both are OS-positive measures satisfying every RP consequence. **So (2.1) cannot
fix the sign of `f₁²`** — and *no* bound built from `dμ_O≥0` alone can. This is T2 §2.1 (`RP is
upper-bound-shaped`) applied to the amplitude. **[T]**

**Width at `β=2.5` (one-plaquette proxy) `[V]`.** Using the exact one-plaquette moments
(`SU(2)` Bessel `M₀(2.5)=0.1341`; `SU(3)` Weyl `M₀(2.5)=0.0524`) and `a₀=s/(4N²)`, `s=1`:
```
f₁²/a₀ ∈ [0, M₀/a₀] = [0, 2.146]  (SU(2)) ,   [0, 1.887]  (SU(3)) .
```
The window is a **factor ≈2** — tighter than "anywhere in `[0,∞)`" — but contains both `0` and `≈1`, so it
is **vacuous for survival**. (The one-plaquette `M₀` is a *proxy* for the full 4D `M₀=Σ_z C(z)`, which is
larger; `[V]`, not the rigorous `[T]` bound (2.1), which uses the full 4D `M₀`.)

---

## 3 · The reduction `[R] ⟺ W` (modulo the normaliser) `[T]`/`[C]`

This is the central adjudication. `W` (`W-WAVE-CHARTER`): `w₁=f₁²/M₀ ≥ w₀>0`. `[R]`: an amplitude floor
`f₁² ≥ δ a₀`, `δ>0`. Since `M₀>0` for a non-degenerate source,
```
(a)  W  ⇒  f₁² > 0            :  f₁² = w₁ M₀ ≥ w₀ M₀ > 0 .                                        (T)
(b)  [R] ⇒ w₁ ≥ δ a₀/M₀ > 0   :  uniform  iff  M₀ ≤ K a₀  (variance bounded above).               (C)
(c)  w₁ = (f₁²/a₀)·(a₀/M₀)    :  a₀/M₀ exact and >0 (one-plaquette).  SU(3), β=2.5:  a₀/M₀=0.530.  (T)
```
So `f₁²>0` and `w₁>0` are the **same qualitative statement** (`sign(w₁)=sign(f₁²)`), and the *uniform*
floors are equivalent exactly when the normaliser is bounded above and below by constants — `a₀ ≲ M₀ ≲ Ka₀`.
But `M₀ = Var_β(O)` is itself a **positivity datum**, and W1's whole point was that **the normaliser is not
the obstruction** (it is positive by construction at every `β`). Therefore `[R]` inherits *precisely* `W`'s
obstruction and **nothing more**: bounding the amplitude at `β=2.5` is bounding the weight at `β=2.5`,
which is the nonperturbative binding/existence-of-atom input. Witness `[4]`. **[T]`/`[C]**

> **What this rules out.** The hope that `[R]` is a *softer*, technical target than `W` (a mere
> "normalization constant" one could bound separately) is **false**: `[R]` and `W` are the same object
> separated by the exact normaliser. There is no shortcut through the amplitude.

---

## 4 · Why `[R]` cannot be bounded at `β=2.5`: five routes, each `[NG]` for its class

We list the five classes of bound that currently exist in the corpus and the exact reason each fails *at
`β=2.5`* (not at strong coupling, where several work). Witness `[2]`,`[3]`.

### 4.1 Kotecký–Preiss absolute convergence `[T]`/`[NG]`
The only **proven** bound on `|c|` is `|log D_fin|≤C·u` with `C` the Ursell constant, finite **iff
`u·K<1`** (`K(β,N)=80e^{2+β}N`, PROOF-P3 §3). At `2.5`:
```
SU(2): c_f(2.5)=2.042, u=1.021, K=1.44×10⁴, uK=1.47×10⁴ ;   KP disk closes at β_KP=0.0067
SU(3): c_f(2.5)=10.64, u=3.546, K=2.16×10⁴, uK=7.66×10⁴ ;   KP disk closes at β_KP=0.0100
```
`β=2.5` is `250–370×` beyond the disk. The decoration gas is **dense** (trivial-character fraction
`0.98→0.17`, `T1`), so this is physical, not a radius-of-convergence accident. **No bound on `c`.**

### 4.2 Reflection positivity / KL positivity `[T]`/`[NG]`
RP gives `dμ_O≥0` ⇒ `C(t)` completely monotone, `log C` convex, `m_eff` decreasing to `m₁` (T2 §2.1,
`BS-BLIND`). These are **upper-bound-shaped**: a pure two-tube **cut** (`f₁²=0`) satisfies *all* of them.
Witness `[3]`: a synthetic atom measure (`f₁²=0.3`, `m₁=0.8`) and a pure-cut measure are **both
log-convex** (`min Δ²log C>0`) with **decreasing** `m_eff` (`→0.800=m₁` vs `→3.366=cut edge`). RP cannot
separate `f₁²>0` from `f₁²=0`. **Sign-blind.**

### 4.3 Exact one-plaquette moments `[T]`/`[V]`
The `SU(2)` (Bessel) and `SU(3)` (Weyl) one-plaquette moments are **exact at every `β`** — but they fix
only `M₀(β)` (the normaliser), i.e. the *denominator* of the pin. They say nothing about the numerator
`f₁²` beyond (2.1). The `β=2.5` window `[0, 1.89a₀]` (§2) **contains both endpoints**. **Sign-blind.**

### 4.4 No FKG / Griffiths order at `β>0`, `N≥3` `[NG]`
A one-sided (Griffiths/FKG) inequality that would push `f₁²` up needs positive plaquette-plaquette
couplings, i.e. **all** character coefficients `c_r(β)≥0`. For `SU(2)` this holds (modified-Bessel ratios
>0). For **`SU(N≥3)` it is not established** (P3-DRESSED §4.1 explicitly flags the "audit-flagged,
not-established for `N≥3` blanket `c_r≥0`"); the `d_{abc}` cubic term already makes the `SU(3)` normaliser
non-monotone (`M₀(2.5)=0.0524 < M₀(0)=0.0556`, W1 §2.2). **No monotonicity tool applies.**

### 4.5 Perron–Frobenius / transfer matrix `[T]`/`[NG]`
PF (Krein–Rutman on the positivity-improving gauge transfer matrix) gives a strictly positive vacuum and
bounds the **diagonal** `⟨Ω|O|Ω⟩`. The connected source `φ=QO|Ω⟩` is **sign-indefinite**, and the overlap
`⟨G|φ⟩` is a difference of two positive-vector flows PF cannot separate from zero (T2 §2.2, `INEQ-R1`).
The transfer matrix gives only a **decreasing upper bound** on `ŵ₁(T)` (T2 §2.3), not a floor. The free
Maxwell field is PF with `f₁²=0`. **Sign-blind.**

> **Synthesis.** Every available bound is one of: (i) valid only inside the KP disk (`β≲0.01`), or (ii)
> one-sided in a direction that **does not reach the numerator `f₁²`**, or (iii) unavailable because the
> necessary sign/monotonicity order (`c_r≥0`, FKG) is absent for `SU(N≥3)` at `β>0`. **None bounds the
> sign of `c` at `β=2.5`.** This is `[NG]` for the class of expansion/positivity/one-plaquette/PF
> functionals, and `[O]` for the irreducible statement that `[R]` = the binding input.

---

## 5 · Acid-test battery (charter §2), all five + continuum, on the amplitude object

| # | acid test | amplitude object `f₁²=a₀(1−cu+…)` | label |
|---|---|---|---|
| 1 | **free massless** `f₁²=0` | no `f⊗f̄→1` binding recoupling ⇒ `a₀=0` ⇒ `f₁²=0`; correctly gets **no** floor | **[T]** |
| 2 | **2D YM** `tr F²=H`, `f₁²=0` | `⟨χ_{1/2}|χ_{3/2}⟩=1.4×10⁻¹⁷`; `f₁²=C₂(R)⟨χ_0|χ_R⟩=0` | **[T]**/**[V]** |
| 3 | **one plaquette** degeneracy | no transfer direction; `a₀/M₀=s/2` is an **artefact** (atom degenerate) — flagged, not invented | **[T]** |
| 4 | **strong coupling** `w₁=s/2` | `f₁²=a₀`, `M₀=1/(2N²)` ⇒ `s/2` (`N≥3`), `s/4` (`SU(2)`) reproduced | **[T]** (mod `[KL]`) |
| 5 | **large `N`** `w₁=O(N⁰)` | `a₀=1/(4N²)`, `M₀=1/(2N²)` both `O(N⁻²)` ⇒ `w₁=O(N⁰)`, `N`-uniform | **[T]** |
| 6 | **continuum** `γ_{w₁}=0` | `f₁²→Z_O²f₁²`, `M₀→Z_O²M₀` ⇒ `w₁` invariant under `O→Z_O O` | **[T]** |

All six pass. Crucially, tests 1–2 confirm the amplitude tool **returns `0`** where the atom is absent —
so it is non-vacuous (it is the *free-field acid test* that the normaliser-only continuation failed in W1).

---

## 6 · Ledger

| item | statement | label |
|---|---|---|
| amplitude coefficient `c` | `O(1)` end-cut Ursell sum in `D_fin=1−cu+…`; **resummed, not finite-order** | **[O]** |
| pin | `0 ≤ f₁² ≤ M₀` exact at every `β`; both ends attained by exact models | **[T]** |
| pin width at `2.5` (proxy) | `f₁²/a₀ ∈ [0, 2.15]` (`SU(2)`), `[0, 1.89]` (`SU(3)`) — sign-blind | **[V]** |
| reduction | `sign(w₁)=sign(f₁²)`; `[R] ⟺ W` uniform iff `a₀ ≲ M₀ ≲ Ka₀` (normaliser not the obstruction) | **[T]**/**[C]** |
| KP bound on `c` | `abs(c) ≤ C` finite iff `uK<1`; at `2.5` `uK=1.5×10⁴` (`SU(2)`), `7.7×10⁴` (`SU(3)`) — **vacuous** | **[T]**/**[NG]** |
| RP / KL | completely monotone, `log C` convex hold at `f₁²=0` too — upper-bound-shaped | **[NG]** |
| one-plaquette | `M₀(β)` exact; fixes denominator only, not `f₁²`'s sign | **[T]**/**[NG]** |
| FKG / `c_r≥0` | not established for `SU(N≥3)` at `β>0`; no monotonicity push | **[NG]** |
| Perron–Frobenius | diagonal/anchor only; connected source sign-indefinite; free field is PF with `f₁²=0` | **[NG]** |
| **verdict** | **`[R]` not boundable at `β=2.5` by expansion/positivity tools; it IS the binding input `W`** | **[O]** |
| Clay | unaffected; W1B confirms the `[R]` = `W` = `I₂` = `H*` = gap ledger, no shortcut | — |

---

## 7 · What W1B establishes, and what it does not

- **Does:** identify `c` as an intrinsically **resummed** object (no finite-order value); prove the exact
  two-sided pin `0≤f₁²≤M₀` and its **sharpness** (both ends realised by exact models); show the `β=2.5`
  window is a factor ≈2 yet sign-blind; prove the **reduction `[R] ⟺ W`** modulo the (non-obstructing)
  normaliser; and show **five** bound routes each fail at `2.5` — KP radius `uK≈10⁴`, RP/KL upper-bound
  shape, one-plaquette denominator-only, missing FKG order for `SU(N≥3)`, PF anchor/diagonal-only.
- **Does not:** give any new lower bound on `f₁²` at `β=2.5`; disprove that *some future* dynamical
  inequality could bound `c`; or move the wall. It **confirms** the wall is one object (`[R]`=`W`) and
  rules out the hope that the amplitude is a "soft" separate technicality — it is the binding input itself.
  Consistent with `T-WAVE-SYNTHESIS` §5: the next objects remain `W` (uniform weight floor) and `I`
  (isolation `ρ²≤C²`), independent and both open. **No QED. Clay 0%.**

*Refs.* `W1-STRONG-ALIGN.md` (§3 coefficient split, §5 conditional bracket); `unconditional/PROOF-P3-UNCONDITIONAL.md`
(§3 KP `K≤6.1×10²N`, §4 dressed amplitude); `P3-DRESSED-MASS.md` (§2.3 `φ_dec=O(u)`, §3.1–3.2, §4.1 `|log D_fin|≤Cu`,
§4.4 free-field acid test); `P3-FOURWALL-DERIVATION.md` (§2 `a₀>0`, four-wall). `T2-TUBE-OVERLAP.md` (§2.1 RP
`[NG]`, §2.2 PF `[NG]`, §2.3 monotone upper bound); `T-WAVE-SYNTHESIS.md` (§0 six verdicts, §5 next objects).
Witness `research/v3/w7/w1b_p3_amplitude.py` (`w1b_p3_amplitude.out`); reuse of `w1_strong_coupling.py`
conventions. **No QED. Clay 0%.**
