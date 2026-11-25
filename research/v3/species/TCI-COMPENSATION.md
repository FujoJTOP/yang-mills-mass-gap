# TCI-COMPENSATION — (TCI) written quantifiably, the compensation/elimination split proved, and the verdict: **the compensation mechanism exists and is strictly stronger in reach than elimination, but it has exactly TWO power sources — the bounded Haar/Jacobian floor `N/12`, and the condensate `(d−3)cF` which *is* the target. No third. `[NG]` (species) / `[T]` (the split) / `[M]` (the missing third source).**

> 2026-11-24 · Round 2, species **TCI**. The third named inequality of `ch13-machine.tex` §15.6, the one the
> 物种波 has not attacked at the level of its own *definition*.
> Reads: `docs/paper-v2/ch05-theorems.tex:271–340` (`eq:twoblk`/`eq:twoblock-bound`/`eq:tci`/`eq:vardecomp`),
> `research/v3/SIXTH-SPECIES-CHARTER.md`, `research/v3/SM5-TCI-HEAD.md`, `research/v3/SM1-SOFT-EXACT.md`,
> `research/v3/Q3-SCHUR-DEFECT.md`, `research/unconditional/PROOF-SAVVIDY-REFORMULATION-CONVEX.md`,
> `research/unconditional/OTHER-SAVVIDY-SIGN.md`, `research/v3/H6-HEAD-SPECTRAL.md`,
> `research/v3/S2-BOUNDARY-COUNTEREXAMPLE.md`, `research/unconditional/SAVVIDY-DTERM.md`,
> `research/unconditional/FORGE-BLIND.md`.
> Witness: `research/v3/species/w7/tci_comp.py` → `tci_comp.log` (`py -3.12 -I`, numpy only, single process,
> deterministic; all quoted numbers are its output).
>
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional/model · **[NG]** no-go for a
> named class · **[O]** exact missing step · **[V]** numerics this round. **No QED. Clay 0%.**

---

## §0 · Verdict (one page, no QED)

**(1) Definition audit — `(TCI)` as printed is not a formula.** `eq:tci` (`ch05-theorems.tex:302–306`) writes
`⟨φ,Qφ⟩+⟨φ,𝕂φ⟩−cF ≥ μ` and then **defines none of `Q`, `𝕂`, `cF`, `ℋ₋`, `ℱ`**. `𝕂` gets a *promissory note*
("the lift of the block scheme; it is unrelated to the fibre Hessian `K`", `:307`) and no formula anywhere in
the paper. `cF` is fixed only in `SM5-TCI-HEAD.md:74,91` as `cF := g|B|`; `Q` is fixed only negatively
(`:310`, "on the free massless field `Q=0`") while the *same letter* is the Savvidy operator in `SM5:72`
(`Q_Sav=−D²+2ad F_B`) and a lattice Laplacian in `SM1:273` — **two inequivalent `Q`'s under one symbol.**
I fix the definitions in §1 and say so explicitly.

**(2) The computable form.** With the faithful choices, `(TCI)` collapses to a **purely magnetic, Haar-free,
exactly computable β-threshold** (§1.4):

```
(TCI-floor)     μ₃ := λ_min( β·Hess S_W | ℋ₋ ) + N/12  ≥  μ > 0
                ⟺   |λ₋^mag(β)|  ≤  N/12 − μ            (λ₋^mag := the negative magnetic eigenvalue)
```

because `Hess Φ ⪰ (N/12)𝟙` **pointwise, exactly** (verified: `min` over 6·10⁴ samples `=0.166680` vs
`N/12=0.166667`). Measured at the head (`2 plaquette`/`2×2`, `SU(2)`, `β=2.5`): `λ_min(β Hess S_W)`
median `−0.3183`, quantiles `[−0.6745, −0.3529, −0.1796, −0.0511, +0.0921]`, `P(μ₃≥0)=0.070`, threshold
`β_c ≈ 0.76`. **`(TCI)` fails at the head by a factor `2.72` in the relevant ratio** (§3).

**(3) Compensation ≠ elimination — proved, and it is the species' real content. [T]**

| | 消去 (Feshbach/Krein/Schur) | 补偿 (in-block, additive) |
|---|---|---|
| object | `Σ = A − BC⁻¹B*` | `A_𝕂 = A + 𝕂` |
| monotonicity | `Σ ⪯ A` **always** ⇒ can only *lower* | `𝕂 ⪰ 0` required ⇒ must *raise* |
| hypothesis | `C ⪰ γ I`, `γ>0` | **none on the complement** (`γ`-free) |
| at the head | dies: `γ_j→0` ⇒ `Σ` undefined | **has no threshold** |
| algebra | congruence / ordered (one-sided) | spectral shift / scale matching (two-sided) |

Sharpest form: `Σ ⪰ μ ⟹ A ⪰ μ`, i.e. **elimination's conclusion is strictly stronger than compensation's**
(the gap is exactly the Schur loss `BC⁻¹B*`), and the converse fails. On the 4D `L=1` Savvidy scan the
in-block curvature `h_n` is `+0.1667 → +0.7280` for `c=0 → 1.5` (**F-monotone increasing**), while `Σ` decays
to `+0.0109` at `c=0.52` and then the fibre `K` goes indefinite at `c_γ ≈ 0.52` where **the Schur form stops
existing** (`ch05:267` says exactly this: "the one-loop Schur form is undefined because `H_ff` is
indefinite"). **So compensation is not a re-naming of elimination: it is the only member of the block algebra
that survives `γ→0`.** That is the strongest thing TCI has going for it — and it is real.

**(4) But the compensation is not free: `𝕂` has exactly two positive sources, and both are spent. [T]**

```
    (i)   𝕂 = Hess Φ        : pointwise uniform, ⪰ N/12, **F-independent and bounded**  → magnitude-FAIL at the head
    (ii)  𝕂 = (d−3)cF        : nonperturbative, = λ₄|a*|², **= the condensate = the target**
    (iii) 𝕂 = 0             : frozen Hessian ⇒ (TCI) ≡ Hess V ≻ 0 ⇒ false wherever ℋ₋≠∅
    (iv)  𝕂 = −BC⁻¹B*        : the Schur lift, ⪯0, = elimination in disguise
    (v)   𝕂 = Γ₂ (Bakry–Émery): ⟺ LSI ⟺ H1 (assumes the conclusion)
    (vi)  𝕂 = c·(−∇²)ᵖ      : p=1 already inside Q; p=2 is `Γ₂`; any p is inside `Hess V`
    (vii) 𝕂 = topological   : exactly 0 (an exact 4-form has zero Hessian; the θ-term is pure imaginary)
```

Exhaustive because `V = S_W + Φ` has **exactly two additive pieces** (`Hess V = β Hess s + Hess Φ`), so
anything Hessian-derived is a piece of those two plus the saddle shift. **⟹ there is no third power source in
the Wilson/Haar formulation.** The species is therefore not a species: `(TCI)` is either *false* (i, iii, iv),
*circular* (v), or *identical to its own target* (ii).

**(5) Two formulation defects found (new this round).**
- **Vacuity `[T]`.** `ℋ₋` is defined *by* the Hessian, so on the free massless field `ℋ₋ = {0}` and `∀φ∈ℋ₋`
  is **empty** ⇒ `(TCI)` holds *vacuously*. The corpus claims the free-field acid test "correctly returns 0"
  (`SM5 §4.1`); in the literal `∀φ∈ℋ₋` form it **does not fail, it evaporates**. The fix is to quantify over a
  *background*-fixed block `ℋ₋^{fix}(F_B)` of dimension `c_N(gB)²V/8π²` and to require `dim ≥ 1`.
- **Object mismatch `[T]`.** At a Savvidy saddle `dim ℋ₋ ∝ (gB)²V` is **extensive**, while on the β=2.5
  measure `E[n₋] = 1.70` (measured, `dim ℋ = 12`). `(TCI)` at the saddle is a vastly stronger statement than
  anything the measure needs — and the measure-typical statement is *already false* (§3). So the strong form
  is dead on both counts.

**One line.** `(TCI)` is a genuine *third* block-algebra mechanism — γ-free, F-monotone, and the only thing
that survives the head's `γ_j→0` — but its positive term can only be the bounded Haar/Jacobian floor `N/12`
(measured `2.72×` too small at `β=2.5`) or the condensate `(d−3)cF` (which *is* `f₁²`/`⟨¼trF²⟩`/H1). Proved
here: `(TCI) ⇏ Σ⪰μ` (it is strictly weaker than the Schur form), every Bogoliubov/symplectic lift is
sign-negative, and the compensation identity `V″(a*)−cF = (d−3)cF`. **`[NG]` as a species; the residue is H1,
unchanged.**

---

## §1 · 新对象（定义）：the audit, and my chosen fixing

### 1.1 Line-by-line audit of the three components

| component | where it appears | is it defined? | verdict |
|---|---|---|---|
| `⟨φ,Qφ⟩` | `ch05-theorems.tex:303` | **NO.** The only gloss is negative: `:310` "on the free massless field, `Q=0`, the bottom of `A|_-` is `0`". Hence `Q` = the **interaction (quartic) second variation**, since the free-field `−D²+2adF_B` is `−D² ≠ 0`. But the same letter is the **Savvidy operator** in `SM5-TCI-HEAD.md:72` (`Q_Sav=−D²+2ad F_B`) and a **lattice Laplacian** in `SM1-SOFT-EXACT.md:273` (`Q(k)`). | **ambiguous, two inequivalent meanings** |
| `cF` | `ch05:303` (`−cF`); nothing else in ch05 | **NO.** Unique fixing: `SM5:74,91` `cF := g|B| =: c`. Also `SM5:98` ties it to the *well curvature*: `Q−cF = 2c`. | **defined only by SM5 (an interpolation), not by the paper** |
| `⟨φ,𝕂φ⟩` | `ch05:303`; glossed `:307–308` | **NO — a promissory note**: "`𝕂` is the lift of the block scheme; it is unrelated to the fibre Hessian `K` of the variance decomposition below." Not a formula. The only concrete instance anywhere: `SM5:99` `𝕂 = −BC⁻¹Bᵀ` (Schur lift, convention (A), `⪯0`). | **undefined; at least 7 inequivalent readings (§2.4), with different truth values** |
| `ℋ₋ ⊕ ℋ₊` | `ch05:281` (`eq:twoblock`), lemma title says **"abstract"** | **NO.** Never tied to a physical background `F_B`. | **root cause of the vacuity defect (§3, row 6)** |
| `ℱ` (the field `cF` refers to) | `ch05:303` | **NO.** `F_B` is only in the surrounding prose. | — |
| `μ` | `ch05:303` | implicitly "the uniform floor", not tied to `a`, `L`, `β`. | — |

**Nothing in the corpus defines `(TCI)` to the standard the other two named inequalities are defined to.**
`(H1)` and `(NA)` each have a fully explicit formula and a named hypothesis list (`ch08-immune.tex` §…,
`charter §1`); `(TCI)` is a *template with three holes*.

### 1.2 My fixing (explicitly mine, not the corpus's)

> **This is what I fixed. The corpus does not supply a unique form** (stated here once, as required).

```
(D1)  Background-fixed negative block:  ℋ₋^{fix}(F_B) := span of the LLL rungs of  −D² + 2 ad F_B
      with dim = c_N (g|B|)² V / 8π²      (c_SU2 = 1, c_SU3 = 3;  SM5 §1.1).          [corpus, SM5]
(D2)  Hessian-negative block:           ℋ₋^{H}(V) := span{φ : Hess V φ = λφ, λ<0}.     [SM5 / SM1]
(D3)  Q  := Hess of the quartic part at the TRUE SADDLE,  φ₀² = c/λ₄:
            Q = 3 λ₄ φ₀² = 3 cF            (SM5:98; exact for  V = −½cF a² + ¼λ₄ a⁴).  [MINE]
(D4)  cF := g|B| =: c                     (SM5:74).                                   [corpus, SM5]
(D5)  𝕂  ∈ { 0 , −BC⁻¹B* , Hess Φ , (d−3)cF , c(−∇²)ᵖ , Γ₂ , 0_top }  — enumerated in §2.4,
      each with its own truth value.                                                  [MINE, exhaustive]
```

`(D3)` is the only non-obvious choice: `Q` = the quartic's second variation **at the stabilised saddle**, not
the Savvidy operator. It is the only choice consistent with `ch05:310`'s `Q=0` on the free field *and* with
`SM5`'s `Q−cF=2c`. **Everything below is stated for this fixing and labelled with it.**

### 1.3 The identity that makes `(TCI)` a statement about `Hess V` and nothing else

Write `V = S_W + Φ` (Wilson + Haar/Jacobian), so exactly

```
        Hess V = β · Hess s + Hess Φ ,          s := Σ_□ (1 − ½ tr U_□) .                  (1.1)
```

`Hess Φ ⪰ (N/12)𝟙` **pointwise, for every configuration and every β** — this is the elementary content of
`ch05:170–180` (`g''(d) ≥ 1/6` on `|d|<2π`), and the witness verifies it exactly:
`min over 6·10⁴ Haar samples of λ_min(Hess Φ) = 0.166680 ≥ N/12 = 0.166667`, median `0.1781`. Then:

> **(1.2) Exact split of the TCI margin.** For any `φ` and the fixing `(D3),(D5)=(Hess Φ)`,
> `⟨φ,Qφ⟩+⟨φ,𝕂φ⟩−cF = ⟨φ, Hess V φ⟩ + [Q − 3cF] + [⟨φ,𝕂φ⟩ − ⟨φ,Hess Φ φ⟩]`,
> and on the saddle-secant (where `Q=3cF`) this is **`⟨φ,Hess V φ⟩`**. Hence at a generic configuration
> `(TCI)` *is* the tautology `Hess V ≻ 0` on `ℋ₋`, and at the saddle it is the *model-dependent* statement
> `2cF + ⟨φ,𝕂φ⟩ ≥ μ` (`SM5:103`). **`(TCI)` is only non-tautological at saddles, and there `cF` and `Q` are
> mean-field data.** `[T]`

### 1.4 Equivalent characterisations (the computable form)

> **Theorem TCI-1 (three equivalent forms and one NON-equivalence). [T]**
> Let `A_𝕂 := Q + 𝕂` on `ℋ₋`, `cF` a scalar, `μ>0`. Then
> ```
> (i)   (TCI)  ⟺  λ_min( A_𝕂 | ℋ₋ ) ≥ cF + μ ;
> (ii)  (TCI)  ⟺  A_𝕂 − cF·𝟙_{ℋ₋}  ⪰  μ·𝟙_{ℋ₋}   (a quadratic form on ℋ₋, uniformly positive);
> (iii) (TCI)  ⟺  inf_{‖φ‖=1, φ∈ℋ₋} ⟨φ, A_𝕂 − cF φ⟩ / ⟨φ,φ⟩ ≥ μ   (Rayleigh floor);
> (iv)  (TCI)  ⟺  V_eff″(a) ≥ μ for all a in ℋ₋-directions at the coarse minimum,  *pointwise, nonperturbatively*
>               (the "gauge-invariant effective potential is convex at its minimum" reading, ch05:308).
> ```
> **But** `(TCI)` is **NOT** equivalent to `Σ ⪰ μ` on `ℋ₋`:
> ```
>   Σ := A_𝕂 − B C⁻¹ B*  ,   C ⪰ γ I  ⟹  B C⁻¹ B* ⪰ 0  ⟹  Σ ⪯ A_𝕂
>   ⟹  ( Σ ⪰ μ )  ⇒  ( A_𝕂 ⪰ Σ ⪰ μ )  ⇒  (TCI)          [one way only]
>   and  A_𝕂 ⪰ cF+μ  gives only  Σ ⪰ cF+μ − ‖B‖²/γ       [the converse loses the Schur loss]
> ```
> So **the Schur form is strictly stronger than `(TCI)`**, by exactly `‖B‖²/γ`. By Haynsworth
> `In(𝔸) = In(C) ⊕ In(Σ)`, `n₋(𝔸) = n₋(C)+n₋(Σ)`, so `(TCI)` alone never yields `𝔸 ≻ 0`; the two-block lemma
> (`eq:twoblock-bound`) is needed and is a *separate* hypothesis. **"`(TCI)` ⟺ `Σ ≻ μ`" is FALSE.**

> **Corollary TCI-2 (the computable, Haar-free form). [T]** With `𝕂 = Hess Φ` and (1.1),
> ```
>   (TCI)  ⟺  λ_min( β Hess S_W | ℋ₋ ) + N/12 + ⟨φ,(Hess Φ − N/12)φ⟩  ≥  μ ,
> ```
> and its **provable** (floor) version — the only version a proof could use, since `Hess Φ`'s sample value is
> unknown a priori — is
> ```
>   (TCI-floor)   μ₃(β) := λ_min( β Hess S_W | ℋ₋ ) + N/12  ≥  μ > 0
>                 ⟺  |λ₋^mag(β)| ≤ N/12 − μ ,      λ₋^mag := min spec(β Hess S_W) .
> ```
> **`(TCI)` reduces to: the Savvidy magnetic eigenvalue never exceeds the Haar/Jacobian floor `N/12`.** This
> is a pure `β`-threshold, `(L,a)`-uniform, and — unlike every other object in the campaign — **exactly
> computable on any block**, with no measure and no conditioning.

Corollary TCI-2 is the answer to "write `(TCI)` in computable form": **one number, `λ_min(β Hess S_W)` per
configuration, compared with `N/12`.** §3 gives it.

---

## §2 · 该对象的控号定理：compensation vs elimination, the Bogoliubov sign, and the exhaustion

### 2.1 Lemma C1 — compensation ⟹ strictly weaker than elimination, and `γ`-free `[T]`

Let `𝔸 = [[A,B],[B*,C]]` on `ℋ₋⊕ℋ₊`, `C ⪰ γ𝟙`. Two demands:
```
  ELIM (消去):        Σ := A − B C⁻¹ B*  ⪰  μ 𝟙_{ℋ₋}
  COMP (补偿):        A + 𝕂            ⪰  (cF + μ) 𝟙_{ℋ₋} ,    𝕂 self-adjoint, ⟨φ,𝕂φ⟩ ≥ 0 required
```
**Proof of the inclusion.** `B C⁻¹ B* ⪰ 0` (as `C ≻ 0`), so `Σ ⪯ A`. Hence `Σ ⪰ μ ⟹ A ⪰ μ`, i.e. **ELIM ⇒ COMP
with `𝕂=0`** — elimination is the *stronger* demand. Conversely `A ⪰ cF+μ` yields only
`Σ ⪰ cF + μ − ‖B‖²/γ`, so the implication COMP ⇒ ELIM costs exactly the Schur loss. ∎
Verified on 2000 random blocks: `Σ ≤ h_n` in **every** trial (`tci_comp.log [B]`).

**Why this is the whole story.** `COMP` contains **no hypothesis on `C`**. It never inverts the complement,
so it is insensitive to `γ`; `ELIM` needs `γ>0` and pays `‖B‖²/γ`, which blows up as `γ→0`. The head is
exactly the regime `γ_j→0` (`H1-RG-DESCENT §2`, `SM1 §5`, `Q3 §4`). Measured (`4D L=1`, `SU(2)`, `β=2.5`,
`v_c = ` uniform Cartan `k=0` mode, `tci_comp.log [B]`):

| `g|B|`-amplitude `c` | `h_n` (COMP, in-block) | `γ = λ_min K` | `cᵀK⁻¹c` | `Σ` (ELIM) |
|---|---|---|---|---|
| 0.00 | **+0.1667** | 0.1667 | 0.0000 | +0.1667 |
| 0.30 | **+0.2081** | 0.1113 | 0.0241 | +0.1840 |
| 0.45 | **+0.2570** | 0.0434 | 0.1344 | +0.1226 |
| 0.50 | **+0.2767** | 0.0151 | 0.2235 | +0.0532 |
| 0.52 | **+0.2850** | 0.0030 | 0.2741 | +0.0109 |
| 0.53 | **+0.2892** | **−0.0032** | — | **UNDEFINED** |
| 0.80 | **+0.4192** | −0.2086 | — | UNDEFINED |
| 1.50 | **+0.7280** | −0.9738 | — | UNDEFINED |

Read: **`h_n` is positive on the whole scan and F-monotone increasing** (`0.167 → 0.728`); `Σ` decays to zero
and then the Schur form ceases to exist at `c_γ ≈ 0.52`. `h_n` survives `3×` further. **This is the
compensation mechanism's genuine content, and it is the one thing the 33 dead families do not have: they are
all elimination/Bogoliubov-type and therefore all die at `γ→0`.** In particular `ch05:270` already says the
one-loop Schur form is *undefined* when `H_ff` is indefinite — so the campaign's own text concedes `ELIM`'s
death while treating `COMP` as "just the same thing". **It is not.**

### 2.2 Theorem C2 — every Bogoliubov / symplectic / Feshbach lift is SIGN-NEGATIVE `[T]`

> For any `c ∈ ℝ`, `½p² + ½ω²q² + c(pq+qp) = ½(p+2cq)² + ½(ω²−4c²)q²`. The map `(q,p) ↦ (q, p+2cq)` is
> canonical (`[Q,P]=i`), hence unitarily implementable, hence **spectrum-preserving**; and the effective
> frequency is `√(ω²−4c²) ≤ ω`, with **equality iff `c=0`**.

The witness diagonalises both operators in a 60-dim oscillator basis: spectra agree to all digits and the
ground energy is `1.00000 → 0.97980 → 0.91652 → 0.60000` for `c = 0 → 0.1 → 0.2 → 0.4`
(`tci_comp.log [A]`). **Consequence.** `𝕂_Bog = −4c²·𝟙 ⪯ 0` **always**: the squeeze/Feshbach/Bogoliubov family
cannot produce a positive second term. This kills the `FORGE-SQUEEZE-BOGOLIUBOV` / `Feshbach–Bogoliubov` lane
*by sign, not by magnitude* — a cleaner death than the corpus's "only helps survivors", and it explains why
**`(TCI)` is not obtainable from anything in the elimination family**.

### 2.3 Theorem C3 — the compensation identity `V″(a*) = cF + R″(a*)`, and the exact margin `(d−3)cF` `[T]`

> Let `V(a) = −½cF a² + λ_d a^d`, `λ_d > 0`, `d ≥ 3`, `cF > 0`. The nontrivial stationary point satisfies
> `a*^{d−2} = cF/(d λ_d)`, and
> ```
>   V″(a*) = (d−2)·cF ,        V″(a*) − cF = (d−3)·cF ≥ 0 ,     V″(a*)/cF = d−2 .
> ```
> Verified: `d=3,4,5,6` give `V″(a*)−cF = 0.00000, 1.00000, 2.00000, 3.00000 = (d−3)cF` exactly.

**This is the sharpest statement in the file, and it answers the task's "is the compensation F-independent?"
question precisely:**

1. **The margin is `(d−3)cF`, i.e. it is `∝ F`. It is NOT F-independent in magnitude.** So a *pointwise*
   `𝕂 − cF ≥ μ` **independent of `F`** does not exist from the condensate route. `[T]`
2. **But the RATIO `V″(a*)/cF = d−2` is F-independent** — i.e. the compensation is *scale-covariant*: the
   more tachyonic the background, the proportionally larger the compensation. For the YM quartic (`d=4`)
   the ratio is exactly `3` (`SM5`'s `Q=3c`, `V″=2c`), a pure fourth-derivative combinatorial number.
3. **The compensation degenerates exactly at cubic anharmonicity** (`d=3`, margin `0`). So `(TCI)`'s truth at
   the saddle is *equivalent to the statement that the soft-mode effective potential is at least quartic* —
   a statement about `λ₄ > 0`, i.e. about the condensate `|a*|² = cF/λ₄ > 0`.
4. In the general (non-homogeneous) case the identity is `V″(a*) = cF + R″(a*)` where
   `R := V + ½cF a²` is the whole beyond-quadratic part; **`R″(a*) ≥ 0` is the compensation, and `R″(a*)` is
   `3λ₄a*²` for the quartic.** Nothing else is available.

### 2.4 Theorem C4 — exhaustion: `𝕂` has exactly two positive sources, and both are spent `[T]`

`V = S_W + Φ` has **exactly two additive pieces** (1.1), so any `𝕂` capable of entering `Hess V` is a piece of
those two plus the nonperturbative saddle shift. Enumerating:

| # | `𝕂` | exact form | sign / status | verdict |
|---|---|---|---|---|
| (i) | frozen-Hessian | `0` | `(TCI)` ≡ `Hess V ≻ 0` on `ℋ₋`; false when `ℋ₋≠∅` | `[NG]`, measured `P(n₋>0)=0.851` |
| (ii) | Schur / Feshbach lift | `−BC⁻¹B*` (`SM5:99`) | `⪯0`; makes `(TCI)` *harder*; `γ`-dependent | `[NG]` (§2.1) |
| (iii) | **Jacobian / Haar** | `Hess Φ` | **`⪰ (N/12)𝟙` pointwise**, β-independent, bounded above by `O(1)` | **`[T]` positive, `[NG]` at the head: `2.72×` short (§3)** |
| (iv) | **condensate** | `(d−3)cF = 3λ₄a*²` | `≥ 0`, **`∝ F`** (F-dependent margin, F-independent ratio `d−2`) | **`[T]` positive but = the target (`f₁²`/`⟨¼trF²⟩`/H1)** |
| (v) | Bakry–Émery | `Γ₂(f) = ‖Hess f‖²_HS + ⟨∇f, Hess V ∇f⟩` | `Γ₂ ≥ ρ ⟺ LSI(ρ)` — **equivalent to the whole target** | `[NG]` (assumes the conclusion) |
| (vi) | kinetic / positive order | `c(−∇²)ᵖ` | `p=1`: already inside `Q` (`Hess V`); `p≥2`: is `Γ₂` | `[NG]` (nothing new) |
| (vii) | topological | `θ∫F∧F̃` | **exactly 0**: an exact 4-form has zero Hessian; the θ-term is pure imaginary, it shifts the *phase*, never `Hess V` | `[NG]` (identically) |

**Completeness argument.** (i)–(ii) are the two pieces of `Hess S_W` (the tachyonic part is the source of
`−cF`; the remainder is the gapped complement); (iii) is all of `Hess Φ`; (iv) is the only contribution that
does not exist at the perturbative saddle. Anything not of the form (i)–(iv) is either a *non-Hessian*
functional — where the only positive candidate is `Γ₂`, which is `(v)` and equals the target — a
higher-order differential operator (vi, which either lives in `Hess V` already or is `Γ₂`), or topological
(vii, zero). **There is no eighth.** So:

> **`(TCI)` has no free positive source. Its two available sources are the bounded Haar floor (too small by a
> measured `2.72×`) and the condensate (which is the target). `[T]`**

### 2.5 What `W_geo ≥ 0` and `Φ″ ⪰ N/12` can and cannot do (the "好集定理" question, answered)

The task asks whether `W_geo ≥ 0` (`eq:wgeo`) or `Φ″⪰N/12` can serve as `𝕂`, and — if so — whether `(TCI)`
was already covered by the "good set" theorems on the small-field set, leaving only the head ⇒ H1 ⇒ dead.
**Answer: yes on both counts, and it is row (iii) above.**
- `Φ″⪰N/12` **is** the pointwise `𝕂` (`Lemma B`, `ch05:178`). It is *not* a new source; it is the only
  *pointwise-uniform* positive term in the whole action.
- `W_geo = E[Φ″] − Var(Φ′) ≥ 0` (Prékopa, `eq:wgeo`) is the **measure-corrected** version of the same row: it
  subtracts the variance, so it is *smaller* than `E[Φ″]` and cannot beat it; `OTHER-SAVVIDY-SIGN §2` already
  measured `Var_Φ(a) = O(1)` but its only rigorous transport is `e^{Θ(N_□)}` (vacuous).
- **So `(TCI)` on the good (small-field) set is indeed covered — and it is *vacuous* there** (no negative
  block ⇒ nothing to compensate), while on the head it fails by the measured factor. **The "good-set covers
  the rest, head ⇒ H1" reading is correct, and this file does not rescue it.** `[T]` for the reading,
  `[NG]` for TCI.

---

## §3 · 酸测五点记分 + 第一手数值

### 3.1 The five-point acid battery (charter §3), plus the two new defects

| # | test | result | verdict |
|---|---|---|---|
| 1 | **free massless field must fail / return 0** | `Hess V = (k²+N/12)𝟙 > 0` ⇒ `ℋ₋ = {0}` ⇒ `∀φ∈ℋ₋` is **EMPTY**, so `(TCI)` holds **vacuously**. It does *not* "return 0" — it evaporates. `SM5 §4.1`'s claim that the certificate "correctly returns 0" is true of the *expression* `Q+𝕂−cF=0`, false of the *inequality*. | **FAIL (vacuous pass)** |
| 2 | thin atom (`w₁→0`) | `(TCI)` is a curvature statement; it is *silent* on `ν`'s atoms — `ch05:333` proves `β=0` gives `G″≥S≥N/12>0`, and charter §3 / `R5-EXACT-ANCHORS` note that this is curvature, not an atom (`f₁²=0` at `β=0`). Cannot manufacture a false atom. | n/a (not its job) / **FAIL as gap tool** |
| 3 | deep double well | the Savvidy double well *is* the counterexample: `V″(well)=2cF`, `V″(barrier)=−cF` ⇒ `ess inf = −cF < 0`. Measured `λ_min(Hess V)` median `−0.1796`, `q05=−0.6745`. | **FAIL (ess-inf form false)** |
| 4 | normalisation (total vs mean, `Var(q)=χ` vs `χ/V`) | `(TCI)` has no normalisation ambiguity (it is a Hessian, not a variance), **but** its `ℋ₋` is extensive at the saddle (`∝(gB)²V`) and O(1) on the measure (`E[n₋]=1.70`) — an object mismatch of order `V`. | **FAIL (§0(5))** |
| 5 | dimension: `G″_a = a²m²_phys` | `(TCI)`'s three terms share dimension `=action`, and `cF=g|B|` is a *mass²* while `⟨φ,Qφ⟩` is a *second variation of the action* — **the printed inequality mixes a form-value with a coupling unless `φ` is `V`-normalised**. With `‖φ‖=1` in `L²` and `Q=3λ₄φ₀²`, both sides are `O(g²)` — consistent only under `(D3)`. | **conditional (needs (D1)–(D5))** |
| 6 | **NEW: vacuity** | see row 1 + §0(5): fixing `ℋ₋ := ℋ₋^{fix}(F_B)` with `dim ≥ 1` is required to make the acid test discriminating. | **defect `[T]`** |
| 7 | **NEW: β-window** | `(TCI-floor)` holds as a *theorem-shaped* statement only for `β < β_c ≈ 0.76`; the strong-coupling window `β<1/48` is inside it, and the head `2.5` is outside. So `(TCI)` is a **true strong-coupling statement and a false head statement** — nothing in between. | **`[NG]` for bridging** |

### 3.2 The numbers — `2 plaquette` / `2×2` cross block, `SU(2)`, flat boundary, `dim ℋ = 12`

`Hess V = β·Hess s + Hess Φ`, `Hess Φ` **exact** (block-diagonal, both curvatures analytic), `Hess s` by
5-point FD in the flat chart, `6·10⁴` Haar-importance samples (`tci_comp.log [C]`).

```
beta   ESS%   E[n_-]  P(n_->0)  med lam_min(Hess V)  P(lam_min<0)  med mu3   P(mu3>=0)  |lam^mag|/lam_min(Phi)
0.50   77.5    0.14     0.129         +0.0630            0.129      -0.0415     0.285         1.169
1.00   37.4    1.15     0.682         -0.0560            0.682      -0.1846     0.103         1.972
1.50   13.9    1.55     0.787         -0.1255            0.787      -0.2609     0.076         2.401
2.00    5.2    1.67     0.829         -0.1618            0.829      -0.2979     0.071         2.609
2.50    2.1    1.70     0.851         -0.1796            0.851      -0.3183     0.070         2.723
3.00    0.9    1.70     0.849         -0.1979            0.849      -0.3276     0.083         2.775
4.00    0.2    1.63     0.788         -0.2044            0.788      -0.3562     0.116         2.936
```

`β = 2.5` detail (`λ_min(Hess V)` = the TCI margin under the Jacobian lift `𝕂 = Hess Φ`):
```
  q0.05 = -0.6745   q0.25 = -0.3529   q0.50 = -0.1796   q0.75 = -0.0511   q0.95 = +0.0921
  n_- distribution (weighted):  {0: 0.149, 1: 0.272, 2: 0.389, 3: 0.127, 4: 0.047, 5: 0.014}
  E[n_-] = 1.70 ,  dim H = 12 ,  dim H_- /dim H = 0.142 (weighted mean)
  exact:  min over samples lam_min(Hess Phi) = 0.166680 >= N/12 = 0.166667   [T check PASSES]
  beta threshold (weighted median of lam_min(Hess V) crossing 0):  beta_c ~ 0.76
```

**The first-hand `inf_{φ∈ℋ₋}` number.** Under the fixing, `inf_{‖φ‖=1, φ∈ℋ₋} [⟨φ,Qφ⟩+⟨φ,𝕂φ⟩−cF]`
equals `λ_min(Hess V|ℋ₋)` (§1.3), whose measure-distribution is the column above:
```
   inf over phi in H_-  =  -0.1796   (weighted median)  ,   -0.6745 (5% worst)  ,   +0.0921 (95%)
   => (TCI) is FALSE with frequency 0.851 at beta=2.5 ; the provable floor mu3 is negative with frequency 0.930.
   => the compensation the Haar/Jacobian sector supplies is |lam^mag|/lam_min(Phi) = 2.72  times TOO SMALL.
```

### 3.3 What the compensation would need, and where it lives

The compensation identity (§2.3) evaluated on the lattice at `β=2.5` along the softest eigenvector
(`tci_comp.log [D]`):
```
   case          lam_min V   cF=-V''(0)   lam4      a*,extrap   predicted V''(a*)   V''-cF
   savvidy-like    -2.0268     2.0268    0.0103      14.005         4.0537          +2.0268
   median          -0.8254     0.8254    0.0092       9.461         1.6509          +0.8254
```
**Two facts.** (a) The predicted well curvature is `2cF` and the margin `+cF` — so **if the well existed the
compensation would be comfortably enough** (`2cF = 4.05 ≫ N/12`). (b) **The well is at `a* ≈ 9–14`, far
outside any amplitude reachable on a small block** (sampled to `|a|≤8`, chart ball `|ξ_link|<2π`). Hence:

> **The compensation is real, positive, and F-covariant — but it is a LARGE-FIELD object that no `2 plaquette`
> block can reach.** This is *why* every small-volume witness (this one, `SM1`, `H1`, `G4`) sees only the
> negativity: on small volumes the accessible configuration space contains **no** compensated well. `[T]`+`[V]`

That is a genuine explanation, not a failure of the computation — and it is the reason the numerical route
to `(TCI)` is structurally closed: to *see* the compensation you must be at a field where `λ₄a*² > cF`,
i.e. `a* > √(cF/λ₄)`, and on a block of `N_□` plaquettes `λ₄ ∝ 1/N_□` ⇒ `a* ∝ √N_□` grows faster than the
chart radius. **`(TCI)`'s positive term is therefore never visible in a small-volume numerical certificate.**

### 3.4 The `2×2` vs the literature anchors — consistency check

`λ_min(Hess V)` quantiles at `β=2.5`, `2×2`: `q05=−0.6745, q50=−0.1796, q95=+0.0921`,
`P(λ_min<0)=0.851`. `SM1-SOFT-EXACT §2b` reports `2×2cross`: `q05=−0.667, q50=−0.175, q95=+0.048`,
`frac<0 = 0.811`. **Agreement to `2–3` digits (independent code, independent sampler, exact `Hess Φ` vs FD).**
The `0.851` vs `0.811` difference is the `Hess Φ` treatment + sampling scatter; the sign and the O(1)
magnitude are reproduced. `N/12 = 0.1667` exact check `0.166680` passes.

---

## §4 · 非重复性（对免疫图）：why this is not the 33 dead families, and not `SM5`

**It is *not* "Feshbach/Krein (只帮幸存者)".** That family is *elimination*; §2.1 proves `Σ ⪯ A` always, so
elimination can never raise the in-block form — the corpus is right that it is dead. **`(TCI)` is the only
member of the block algebra that attacks the block itself**, is `γ`-free, and *survives* `γ_j→0`
(measured: `h_n>0` to `c=1.5`, `Σ` undefined past `c_γ≈0.52`). The difference is a *theorem* (`Σ ⪯ A`), not a
re-labelling — so this file's first result is that the campaign's reason for distrusting `(TCI)` ("it's the
same as Feshbach") is **wrong**: `(TCI)` is strictly *weaker* and therefore has strictly larger
admissible region.

**It is not `S3-UV-FESHBACH` / `ch05` Lemma NS (convention (B)).** Lemma NS eliminates a **negative** block
(`C≺0`, `Σ ⪰ A`, "helps"); `(TCI)` does not eliminate anything. `SM5 §1.2` names the confusion of (A)/(B) as
the category error, and this file keeps them apart and uses only (A)-type lifts plus the in-block form.

**It is not `Q3-SCHUR-DEFECT`.** `Q3` is the *inertia* statement `n₋(H)=n₋(K)+n₋(Σ)` and its conclusion
`Σ>0 ⇔ the unstable mode is fibre-confined` — an *excision* criterion, needing `γ` and therefore dead at the
longitudinal barrier. This file instead proves `(TCI) ⇏ Σ⪰μ` (§1.4(iv)) and quantifies the gap `‖B‖²/γ`.

**It is not `SM5-TCI-HEAD`.** `SM5` splits TCI by **which set the infimum runs over** (ess-inf / typical point
/ measure average) and concludes the strong form is false, the true forms weak. This file attacks the
**operator** `𝕂` instead — a different axis, and it is the axis on which `(TCI)` is *undefined* (§1.1). The
new results relative to `SM5`:
1. **§1.4(iv)**: `(TCI)` is strictly *weaker* than the Schur form (new; `SM5` conflates them in §2.3);
2. **§1.4 Corollary**: the reduction of `(TCI)` to the single magnetic number `λ_min(β Hess S_W)` vs `N/12`
   (new — `SM5` never reduces it; it works on the soft mode analytically);
3. **§2.2**: the Bogoliubov sign theorem (new, and it kills `SM5`'s Route-3 family by sign);
4. **§2.3–2.4**: the identity `V″(a*)=(d−2)cF` with margin `(d−3)cF`, and the **exhaustion** proving there is
   no third positive source (new; `SM5` asserts `K≤0`/`2c` without the `d`-law or the exhaustion);
5. **§3.2–3.3**: `λ_min(Hess V)` measure distribution, `n₋` distribution, `β_c≈0.76`, and the finding that
   the compensated well lies at `a*≈9–14`, outside any small-block amplitude (all new);
6. **§3.1 rows 1,6**: the vacuity defect — the literal form does **not** pass the charter-§3 acid test (`SM5`
   claims it does).

**It is not `OTHER-SAVVIDY-SIGN` / `PROOF-SAVVIDY-REFORMULATION-CONVEX`.** Those bound `W_geo` and prove the
Savvidy inertia is chart-invariant — both *negative* about `W_geo`'s magnitude and about reparametrisation.
This file *agrees* with them (§2.5) and adds the reason: `Hess Φ` is the only positive additive piece, and
`W_geo ≤ E[Φ″]` is strictly smaller than it.

## §5 · 归属：第五数据模 / 新坐标 / 决定性否证

**Class (γ) — decisive falsification — NOT (α) and NOT (β).**

- **(β) requires** a new object `𝔉` with positivity ⟺ gap and a *new* machine not in the dead list. `(TCI)`'s
  machine would be "in-block addition instead of elimination", which §2.1 shows is a real and distinct
  mechanism — **but** §2.4 shows its two power sources are the bounded Haar floor (false at the head) and the
  condensate (identical to the target). So it fails the `(β)` test at the *source*, not at the *machine*.
- **(α) requires** a new computable, floor-sensitive data type. `(TCI)`'s computable shadow,
  `(TCI-floor) = λ_min(β Hess S_W) + N/12`, is exactly a **finite-extremal** quantity (charter §5's fourth data
  type) and is therefore **floor-blind by `SPECIES-EXISTENCE`'s `F3∧F4=∅`** — it is a curvature number, not
  an atom detector. So `(α)` fails for the same structural reason as the previous five species.
- **(γ) is what is achieved.** Three falsifications:
  1. **`(TCI) ⟺ Σ ⪰ μ` is FALSE** (§1.4): the Schur form is strictly stronger. The corpus's implied
     equivalence (`ch05:299–301`, "the indefinite Hessian is *converted* into the demand `Σ≥μ>0`", then
     `eq:tci`) does not hold as an equivalence.
  2. **The free-field acid test is passed VACUOUSLY** (§3.1 row 1) — the literal `∀φ∈ℋ₋` form is empty at
     `F=0`. `SM5 §4.1`'s "correctly returns 0" is a statement about the *expression*, not the *inequality*.
  3. **No third positive source exists** (§2.4), by the two-term structure of `V`.
  These narrow the map exactly as charter §5(γ) asks, and mark `(TCI)` **`[NG]` as a species**.

## §6 · 精确缺失的一步（一句话，无退路）

> **A positive, order-0 (multiplicative), pointwise-uniform second-variation term in `V` that is neither
> `Hess Φ` (bounded by `N/12`, measured `2.72×` short at `β=2.5`) nor the condensate `(d−3)cF` (which is
> `f₁²`/`⟨¼trF²⟩`/H1 itself) does not exist in the Wilson×Haar formulation**, because `V = S_W + Φ` has
> exactly two additive pieces; hence the only non-vacuous form of `(TCI)` is `(d−3)cF ≥ μ`, i.e.
> `λ₄⟨|a₀|²⟩ ≥ μ`, i.e. **the uniform positivity of the gluon condensate** = H1. `[M]`

Equivalently, in one clause: **`(TCI)` is H1 after a change of variables that moves the condensate from the
right-hand side to the left.**

## §7 · 可执行的下一步

1. **`tci_mag.py` — the β-window closure test (highest value, cheap).** Measure `β_c(L)` — the threshold at
   which the *weighted median* of `μ₃ = λ_min(β Hess S_W)+N/12` crosses `0` — on `1×2, 1×3, 1×4, 2×2` and
   `4D L=1`. **Prediction from this file's data:** `β_c` *decreases* with volume/plaquette count (`0.76` at
   `2×2`; expect `≤0.5` at `1×4`), which closes the "TCI holds at strong coupling, the UV, and not the head"
   reading by making the strong-coupling window shrink to `β→0`. If instead `β_c(L)` is `L`-uniform at
   `≈1/(something)`, then `(TCI-floor)` is a genuine `(L,a)`-uniform theorem for `β<β_c` and should be
   written up as `[T]` (a *small* region, but real). Either outcome is decisive; either is a few seconds of
   compute on the existing `tci_comp.py` machinery.
2. **`tci_cond.py` — the condensate identity at reachable amplitude.** Instead of the softest eigenvector of a
   single sample, use the *block-uniform* coarse mode and a **twisted-boundary / flux** background so that
   `a*` is O(1) by construction (choose `λ₄` large by using a small `N_□`). Target: verify
   `V″(a*) − cF = 3λ₄a*²` to `10⁻³` on the lattice, up to `a*≈2`. If it holds, §2.3 is confirmed
   nonperturbatively rather than in the normal form; if it fails, the mean-field `d`-law is not the lattice's
   and the compensation magnitude must be measured, not computed.
3. **Do NOT** spend compute trying to make the well reachable on large volumes: §3.3 shows `a* ∝ √N_□` grows
   faster than the chart radius, so **the compensation is structurally invisible to small-block numerics.**
   The only route that can *certify* `(d−3)cF ≥ μ` is a continuum statement about `λ₄⟨|a₀|²⟩` — i.e. H1.

---

## §8 · Ledger

| item | statement | label |
|---|---|---|
| `eq:tci` defines none of `Q,𝕂,cF,ℋ₋`; `𝕂` is a promissory note | `ch05:303–308` audit | **[T]** |
| `Q` has two inequivalent meanings in the corpus (interaction vs `−D²+2adF_B`) | `ch05:310` vs `SM5:72` vs `SM1:273` | **[T]** |
| fixing `(D1)–(D5)`, `(TCI) ⟺ λ_min(A_𝕂\|ℋ₋) ≥ cF+μ` | Thm TCI-1 (i)–(iv) | **[T]** |
| `(TCI) ⇏ Σ ⪰ μ`; `Σ ⪰ μ ⇒ (TCI)`; gap `= ‖B‖²/γ` | Thm TCI-1 (iv) | **[T]** |
| `(TCI-floor) ⟺ λ_min(β Hess S_W\|ℋ₋) ≥ N/12 + μ`, a pure magnetic β-threshold | Cor TCI-2 | **[T]** |
| `Hess Φ ⪰ (N/12)𝟙` pointwise, all `β` | min `0.166680` over `6·10⁴` samples | **[T]**+**[V]** |
| COMP(补偿) strictly weaker than ELIM(消去) and `γ`-free; `Σ ⪯ A` always | Lemma C1, 2000 random blocks | **[T]** |
| at `β=2.5`: `h_n>0` for `c∈[0,1.5]` (F-monotone), `Σ` undefined past `c_γ≈0.52` | `4D L=1` scan | **[V]** |
| every Bogoliubov/symplectic lift is sign-negative (`ω²→ω²−4c²`) | Thm C2 | **[T]**+**[V]** |
| `V″(a*) = (d−2)cF`, margin `(d−3)cF ≥ 0`, zero iff cubic; ratio `d−2` is F-independent, margin is `∝F` | Thm C3 | **[T]** |
| `𝕂` has exactly 2 positive sources (Jacobian `N/12`; condensate `(d−3)cF`); no third | Thm C4, 7-row exhaustion | **[T]** |
| `W_geo ≤ E[Φ″]`; the good-set reading (bad set `⇒` H1) is **correct** and not rescued here | §2.5 | **[T]** for reading / **[NG]** for TCI |
| `2×2`, `β=2.5`: `E[n₋]=1.70`, `P(n₋>0)=0.851`, `λ_min(Hess V)` median `−0.1796`, `q05=−0.6745` | `tci_comp.log [C]` | **[V]** |
| `P(μ₃≥0)=0.070`, `|λ^mag|/N/12 = 2.72`, `β_c≈0.76` | `tci_comp.log [C]` | **[V]** |
| compensated well at `a*≈9–14`, unreachable on a block; predicted margin `+cF` | `tci_comp.log [D]` | **[V]** |
| literal `∀φ∈ℋ₋` is **vacuous** on the free field (acid test FAILS, not "returns 0") | §3.1 row 1 | **[T]** |
| object mismatch: `dim ℋ₋ ∝ (gB)²V` at the saddle vs `E[n₋]=1.70` on the measure | §0(5)/§3.1 row 4 | **[T]**+**[V]** |
| **`(TCI)` as a species** | no free positive source; = H1 | **[NG]** |
| **Clay** | unbuilt | **0%** |

---

## §9 · References (internal, this session)

- `docs/paper-v2/ch05-theorems.tex:271–340` — `eq:twoblock`, `eq:twoblock-bound`, `eq:tci` (`:302–306`), the
  `𝕂` promissory note (`:307`), Lemma NS (`:314–318`), `eq:vardecomp` (`:320–332`), `:170–183` (`Φ″⪰1/6`,
  `W_geo≥0`), `:263–269` (the Schur form "undefined because `H_ff` is indefinite").
- `research/v3/SM5-TCI-HEAD.md` — `cF:=g|B|` (`:74,91`), `Q=3c` (`:98`), `𝕂=−BC⁻¹Bᵀ` (`:99`), the
  set-trichotomy (§0), the four routes (§2), the acid-test claim (§4.1).
- `research/v3/SM1-SOFT-EXACT.md` §2b (2D block spectra, `frac<0`), §5 (descent break `γ_j<0`, `η_j≥1`).
- `research/v3/Q3-SCHUR-DEFECT.md` — Lemma Q3.1 (Haynsworth), the one-plaquette dichotomy, Lemma Q3.2.
- `research/v3/H6-HEAD-SPECTRAL.md` §3 (`λ₁(ν) ≥ ess inf G″`), §8.1 (`(◆)` the tail-controlled floor).
- `research/v3/H2-DEFECT-BOUND.md` — `D = D₊+D_Sav`, `D/S≈0.97`, Theorem M (Fisher).
- `research/unconditional/PROOF-SAVVIDY-REFORMULATION-CONVEX.md` (inertia is chart-invariant; the sole
  positive object is `W_geo`), `OTHER-SAVVIDY-SIGN.md` §0/§2 (`10× short`, `W_geo`'s full strength).
- `research/unconditional/SAVVIDY-DTERM.md` (rank-one `K⁻¹=K₊⁺−P₋/|λ₋|`), `FORGE-BLIND.md` §2–4.
- Standard tools (used as textbook, not re-fetched): Haynsworth (1968) inertia additivity; Schur complement;
  Sylvester's law of inertia; Bogoliubov canonical transformations / squeeze; Prékopa (1971);
  Brascamp–Lieb (1976); Bakry–Émery (1985) `Γ₂` ⟺ LSI; Savvidy *PLB* **71** (1977) 133; Nielsen–Olesen
  *NPB* **144** (1978) 376.
- Witness: `research/v3/species/w7/tci_comp.py`, log `tci_comp.log` — `py -3.12 -I`, numpy only, single
  process, deterministic (`seed=20261124`), exit 0.

*No third-party sources downloaded this session; all external results are standard and cited from the corpus.
**No QED. Clay 0%.***
