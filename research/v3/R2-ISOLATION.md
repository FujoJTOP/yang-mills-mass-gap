# R2-ISOLATION — attacking the atom's missing half `I2`: isolation is **not** a consequence of the atom, it is the **primitive** the atom derives from; moments/positivity are **blind** to it (thin-gap `[NG]`); it is exact in every solvable confining model and is the **order-2 (second-eigenvalue) rung** above the atom's order-1 (Perron)

> **2026-11-24 · Wave R2 — the isolation half of `H*`.** Target: `f₁²=|⟨0|tr F²|G⟩|² ≥ r₀>0`
> **with isolation** `I2: Λ★−m₁²≥δ>0`, `Λ★:=inf supp(ν−f₁²δ_{m₁²})`, `ν` the Källén–Lehmann
> measure of the gauge-invariant `0⁺⁺` composite `tr F²`. Reads `Q-WAVE-SYNTHESIS.md`,
> `Q6-H1-STATUS.md`; `unconditional/PROOF-NO-DISSOLUTION-KINEMATIC.md`, `PROOF-ISOLATION.md`,
> `PROOF-EDGE-E2.md`; and the banked isolation wave `G2D/G2F/G2G/G2H`, `STATE-DENSITY-I2`,
> `SUB2-I2-STRONG`, `SUPPORT-NOT-ATOM`, `INEQ-T1`, `INEQ-T5`; `docs/paper-v2/ch02/ch03/ch05`.
> Witness: `research/v3/w7/r2_isolation.py` (W1–W6, run `py -3.12`, exit 0; all PASS).
> **Labels (charter §2).** **[T]** theorem · **[C]** conditional · **[NG]** no-go for a named
> class · **[O]** exact missing step · **[V]** numerical witness. **No QED; Clay 0%.**

---

## 0 · Verdict (up front)

**Question asked: does isolation `I2` follow from the same structure as the atom `I1`, or is it
genuinely independent?** The wave's answer, sharpened to three statements:

1. **Logically, isolation is the primitive and the atom is its consequence — the reverse of the
   corpus's framing `[T]`.** `I2 ⟹ I1` (a support gap forces `ν({m₁²})>0`; Lemma A, §2); the
   converse `I1 ⟹ I2` is **false** (a continuum may touch `m₁²` from above; Lemma B). So the pair
   `(I1∧I2)` in `(T★)` is, up to one side-condition, **`I2` alone**. The atom is *downstream* of the
   isolation, not the other way. (This corrects Q6 §5's "the irreducible single object is the atom":
   the irreducible object is the **isolation**; the atom is what it hands you.)

2. **At the positivity / moment / sum-rule level, isolation is unreachable `[NG]`.** Moments give
   only an **upper** bound on `Λ★` (Lemma C, §2), exactly mirroring the Christoffel **upper** bound
   on `f₁²` (`INEQ-T1`). The lower end of the admissible isolation is **0**, and it is *realised*:
   the **thin-gap** measures of §3 match the first `K` moments of an isolated spectrum but drive
   `Λ★−m₁²→0` (witness W4: identical `(M₀,M₁,M₂)=(2,11,101)`, isolation `9→0.01`). So `I2` is at
   least as hard as `I1`; no positivity functional lower-bounds it.

3. **Structurally, isolation sits one rung higher than the atom: order-2 vs order-1 `[T]+[NG]`.**
   The atom is an **order-1 positivity** statement (Perron/Krein–Rutman on the single positive
   cone: a strictly-positive transfer kernel has a simple top eigenvalue with a positive
   eigen*vector*, and the overlap `f₁²=N_O(z₁)≠0`). Isolation is a **second-eigenvalue** statement
   `λ₁−λ₂>0` — an **order-2 (total-positivity / second-invariant-cone)** datum. Every order-1 tool
   (PF, RP cone, Bushell, Bakry–Émery, Cheeger) is, by `PROOF-EDGE-E2`, **provably blind to
   `λ₁−λ₂`** (`{1,λ₁,λ₁,…}` and `{1,λ₁,λ₂,…}` share every order-1 datum). This *explains* the
   banked no-gos as a **tier**, not a constant: isolation needs the second Karlin rung, which
   `INVENT-SECOND-CONE`/`H6 §8` show is **absent on the orderless quotient `𝒜/𝒢`**.

**Net.** Isolation is a **genuinely independent input**, strictly stronger than the atom, one
`order-2` tier above it. It is **exact and positive in every solvable confining model** (2D YM:
Casimir ratio `ρ²=8/3`; 1D rotor: exact gaps) and **correctly absent** on the free massless `tr F²`
(no atom, support from 0 — acid PASS). In 4D it is `[O]`, and — via single-scale covariance — the
*fixed-absolute-`δ`* form is **equivalent to the gap itself**, with the *scale-invariant* content
being the dimensionless ratio `ρ=m_{0⁺⁺*}/m₁≥1+c`. **The atom does not give isolation; isolation
gives the atom.**

---

## 1 · Survey of the banked isolation lemmas, and the exact step each names

Isolation has been attacked six times; the wave consolidates them. Every one lands on the same
object in a different coordinate — this is the `I2`-analogue of the Q-wave "four faces" of `H1`.

| route | object | verdict | exact missing step |
|---|---|---|---|
| **No-dissolution** (`PROOF-NO-DISSOLUTION-KINEMATIC`) | kinematic pin `m_{0⁺⁺}<2m_g` | `[C]` | needs the gap + ordering (K-A,K-B); margin `=m_g`, inherits the gap |
| **Edge E2 / Temple** (`PROOF-EDGE-E2`) | second level `Λ=E₂` | `[NG]` | bottom-Dirichlet functionals blind to `λ₁−λ₂`: the Temple regress |
| **Threshold sharpening** (`PROOF-ISOLATION §4.4`) | operative `Λ`=excited `0⁺⁺*≈1.5m_g` | `[T]` (structure) | `C/P` cannot separate `gg` from `(gg)(gg)`; `2m_g` overshoots |
| **SSC / G2D** | ratio `ρ=m_{0⁺⁺*}/m₁` | `[T]`+`[NG]` | fixed-`δ` `I2` ⟺ gap (one-scale); escape = scale-invariant `ρ≥1+c` = same `O(u)` wall |
| **Determinant / G2F** | `D_O(z)` first zero | `[NG]` | zero sits *on* the orbit-radius `z₁`; need `|D_O/D_lead−1|=O(u)` past the edge |
| **State density / STATE-DENSITY** | window count `N(L,a)` | `[T]`+`[NG]` | window count is `O(1)`, `L`-uniform (Weyl benign); residue is `ρ` (`O(u)` binding) |
| **Hankel / INEQ-T1, sum-rule / INEQ-T5** | `f₁²` from moments | `[NG]` | Christoffel is an **upper** bound; literature: "positivity forces minimal residue to zero" |

The four with the sharpest phrasing agree by different routes: **the fixed-`δ` isolation is
equivalent to the gap (SSC), and its scale-invariant content is `ρ≥1+c`, bounded by the same
`O(u)` cost-vs-rate wall that obstructs the atom.** So the wave's task reduces to: *is `I2`
strictly harder than `I1`, or the same object renamed?* §2–§3 answer this cleanly, and give the
first **two-sided logical** separation plus a **model-positive** result.

---

## 2 · The logical lattice: `I2 ⟹ I1`, `I1 ⇏ I2`, moments upper-bound `I2`

Throughout `ν≥0` on `[m₁²,∞)`, `m₁²:=inf supp ν`, atom `f₁²:=ν({m₁²})`,
`Λ★:=inf supp(ν−f₁²δ_{m₁²})`, moments `M_k=∫s^k dν`, mean `R̄=M₁/M₀`.

> **Lemma A (isolation forces the atom) `[T]`.** If `Λ★>m₁²` then `f₁²>0`.
> *Proof.* If `f₁²=0` then `ν−f₁²δ=ν`, so `Λ★=inf supp ν=m₁²`. Contrapositive: `Λ★>m₁² ⇒
> f₁²≠0`, and `f₁²≥0` gives `f₁²>0`. ∎

> **Lemma B (the atom does not force isolation) `[T]`.** There is `ν` with `f₁²>0` and
> `Λ★=m₁²`. *Proof.* `ν=f₁²δ_{m₁²}+μ`, `μ` a continuous density on `[m₁²,∞)` with `μ([m₁²,m₁²+ε])>0`
> for every `ε>0` (e.g. `dμ=1_{(m₁²,∞)}ds`). Then `Λ★=m₁²`. ∎

So **`I2 ⊋ I1`**: isolation is strictly stronger. Consequence, recorded as a correction:
the target `(T★)=(f₁²≥r₀ ∧ I2)` is, qualitatively, **`I2` alone** (`I2 ⟹` atom positive). The
*quantitative* floor `f₁²≥r₀` needs in addition the **weight-dominance** side condition
`Λ★>R̄` (below), because of the heavy-tail escape; but the atom is never an *independent* conjunct.

> **Lemma C (moments upper-bound isolation) `[T]`.** For any `ν` with atom `f₁²` at `m₁²` and
> `M₀>f₁²`,
> ```
>    Λ★  ≤  Λ_max(f₁²,M₀,M₁) := (M₁ − f₁² m₁²)/(M₀ − f₁²).
> ```
> *Proof.* `M₁=f₁²m₁²+∫_{≥Λ★}s dμ ≥ f₁²m₁²+Λ★(M₀−f₁²)`. ∎ (Validated: W6 — for a two-atom
> measure `Λ_max` equals `Λ★` exactly; with a heavy tail `Λ_max` over-shoots, so it is *upper*-only.)

Lemma C is the **exact dual** of `INEQ-T1`'s Christoffel bound `f₁²≤1/K_n(x,x)` on the atom:
*moments bound `f₁²` above and `Λ★` above; neither has a lower companion.* The residue–level–gap
lemma is the same inequality rearranged from the other side (isolation ⟹ atom):

> **Lemma D (residue–level–gap; isolation ⟹ atom floor) `[T]`.** If `Λ★>R̄` then
> ```
>    f₁²  ≥  M₀ (Λ★ − R̄)/(Λ★ − m₁²)  >  0 .
> ```
> *Proof.* Reproduced in `ch05 §sec:reslevel`; equality iff the continuum is one atom at `Λ★`. ∎
> **Side condition `Λ★>R̄` is load-bearing:** a positive atom with a heavy tail can have `R̄>Λ★`,
> in which case `Λ_max<R̄` and the rearrangement gives a *negative* "bound" — the atom is still
> positive (Lemma A) but the *quantitative* floor is lost.

**Logical heart.** The chain is a **one-way diamond** with `I2` the apex-source and `I1` the sink:
```
 I2 (support gap)  ──Lemma A──►  I1 (atom>0)          (qualitative, always)
 I2 + (Λ★>R̄)      ──Lemma D──►  I1 with f₁²≥r₀        (quantitative, needs weight dominance)
 I1                ──Lemma B──X►  I2                    (FALSE: continuum may touch m₁²)
 f₁², M₀,M₁        ──Lemma C──►  Λ★ ≤ Λ_max             (moments UPPER-bound isolation)
```
The atom is **not** an independent primitive; isolation is. This is the clean answer to the
task's question, and it is a **correction** to Q6 §5's "the atom is the irreducible object."

---

## 3 · The new argument, and where it dies: the thin-gap `[NG]` and the finiteness theorem

### 3.1 Thin-gap `[NG]` — moments/positivity cannot lower-bound `I2` (witness W4)

`INEQ-T1` proved the atom is moment-blind (atom weight `0` and `1` share two moments). The **same
obstruction recurs for isolation**, one level up. Take the isolated reference `ν₁=δ₁+δ₁₀`:
`(M₀,M₁,M₂)=(2,11,101)`, isolation `=9`. For any `ε>0` build three atoms `{1,1+ε,L}` with weights
`{p,q,r}` chosen to reproduce `(M₀,M₁,M₂)` **exactly**; such a measure has `Λ★=1+ε`, i.e.
**isolation `=ε`** by construction. Witness W4:

```
   eps    p@1     q@(1+eps)  r@L      L      M1-chk    M2-chk   isolation
  0.500  0.82149  0.20000   0.97851 10.096  11.000000 101.000000   0.5000
  0.050  0.98022  0.02000   0.99978 10.001  11.000000 101.000000   0.0500
  0.010  0.99601  0.00400   0.99999 10.000  11.000000 101.000000   0.0100
```
Identical protected moments, isolation `0.5→0.01→0`. By the Carathéodory count, *any fixed finite
set* of moments is matched by an isolation-`→0` measure (**the continuum can always be pushed up
against `m₁²` and its low moments compensated by a larger mass higher up**). So:

> **`[NG]` (moment/positivity class).** *No functional of finitely many moments, and no
> positivity/Hankel/Stieltjes datum, lower-bounds `Λ★−m₁²`.* The lower end of the admissible
> isolation is `0`, and it is realised. (This is `INEQ-T1` transplanted to the second level;
> `INEQ-T5`'s literature scan — "positivity forces the minimal residue to zero" — is the same wall.)

### 3.2 The finiteness theorem — isolation IS moment-extractable iff the tower is finite (W5)

The thin-gap escape shows moments fail **because the atom count is unbounded**. If the spectrum is
known a priori to be a **finite `K`-atom tower**, the first `2K` moments *determine* it (Prony /
Padé), and isolation becomes an **algebraic function of finitely many moments**:

> **Theorem (Finite-tower isolation) `[T]`.** Let `ν=Σ_{i=1}^K w_iδ_{x_i}` with `K<∞`. Then
> `(x_i,w_i)` are recovered from `(M₀,…,M_{2K−1})` by Prony, and `Λ★−m₁²=x_{(2)}−x_{(1)}` is an
> exact algebraic function of those moments.
> *Proof.* Prony: the atoms are the roots of the degree-`K` polynomial with coefficient vector
> solving the `K×K` Hankel system `Σ_j a_j M_{n+j}=−M_{n+K}`; weights from the Vandermonde. ∎
> **W5** recovers `{1,2.5,7}` from its six moments to `|err|=1e-14`, isolation `1.5` exact. With an
> **8-atom continuum** probed at the *wrong* `K=3`, the recovered "isolation" `1.88` is meaningless.

**Reading.** In solvable confining models the 0⁺⁺ tower below threshold is **finite** (`K` small),
so isolation is moment-accessible and *exactly computable*. In 4D the count of 0⁺⁺ levels in the
window `(m₁²,Λ★)` is **also** finite (`STATE-DENSITY-I2`: an `O(1)` integer, `L`-uniform), **but
`K` itself is the unknown** — precisely the state-density datum `I2` names. So the finiteness
theorem **relocates** the missing input to "know the atom count below threshold", i.e. to `I2`,
rather than removing it. It is the cleanest *positive* statement available: *isolation is algebraic
in the moments conditional on the state count, and the state count is the open datum.*

### 3.3 The structural tier — isolation is order-2, the atom is order-1 `[T]+[NG]`

Why do **every** order-1 positivity tool (PF, RP cone, Bushell/Krein–Rutman, Bakry–Émery, Cheeger)
fail on isolation *for the same reason*? Because they are **order-1 statements about the top of
the spectrum** (a single positive cone, Perron vector, strictly-positive kernel ⇒ simple top
eigenvalue). The atom is exactly such an order-1 datum: `f₁²=N_O(z₁)≠0` with `z₁` the *simple top
zero* of `D_O` — positivity-improving reaches it (`PROOF-GAUSS-INVARIANT-PF`). But isolation is a
**second**-eigenvalue statement `λ₁−λ₂>0`; the invariances of order-1 theory —
`{1,λ₁,λ₁,λ₁,…}` vs `{1,λ₁,λ₂,λ₂,…}` — share every order-1 datum, so **no order-1 cone sees
`λ₁−λ₂`** (`PROOF-EDGE-E2`, `G2G §4`). Seeing `λ₂` requires a **second invariant cone below the
first** — an **order-2 total-positivity / Karlin-flag** structure (a Sturm/oscillatory or
Jacobi/banded transfer), which `INVENT-SECOND-CONE` and `H6 §8` show is **absent on `𝒜/𝒢`**
(the gauge quotient carries no order). 

This is the wave's structural finding: **the atom and the isolation are separated by a Karlin
tier.** Order-1 positivity gives the atom (and the *bottom* gap `E₁−E₀`); order-2 total positivity
would give the isolation (`E₂−E₁`), and it is precisely the second cone the orderless quotient
lacks. This *unifies* the banked no-gos as a tier (not a constant), and predicts they fail
*together* — which is what the SSC equivalence below encodes.

---

## 4 · Model tests (exactly solvable): isolation is real, computable, and acid-clean

The decisive test of any claim is the model battery; all numbers below are from
`research/v3/w7/r2_isolation.py` (W1–W3, deterministic). **The free massless `tr F²` must acquire
no atom; it does not.**

**(W1) Free field, `tr F²` carrier — ACID PASS `[T]+[V]`.**
```
  massless trF^2 (2-gluon cut)      f1^2=0.000  inf supp=0.000   no atom
  massive :phi^2: (cut from (2m)^2) f1^2=0.000  inf supp=4.000   no atom (gapped, atomless)
  elementary free phi (delta_m^2)   f1^2=1.000  inf supp=1.000   atom, but NO continuum
```
The composite free field has `f₁²=0` in both cases and **no isolation** (support is a cut touching
its own bottom; `Λ★=m₁²`). The elementary `φ` is trivially isolated (single atom, no continuum) but
is *not* the `0⁺⁺` composite. **Any instrument that certifies an isolated atom on the free massless
`tr F²` is refuted.** Acid PASS, non-vacuously (it also handles the atomless gapped `:φ²:`).

**(W2) 2D YM (Migdal) — isolation is a symmetry-fixed Casimir ratio `[T]+[V]`.** The 2D transfer
matrix is `e^{−(g²/2)C₂}`, so the 0⁺⁺ tower is `{C₂(R)}`:
```
  SU(2):  C2(j)=j(j+1) for j=1/2,1,3/2,2  ->  3/4, 2, 15/4, 6
  lightest m1^2 ∝ 3/4 ; next Λ^2 ∝ 2 ;  ρ^2 = Λ^2/m1^2 = 8/3 = 2.6667 (EXACT)
  isolation margin Λ^2 − m1^2 = 5/4 > 0  ->  I2 HOLDS exactly
```
In 2D YM the isolation is a **pure group-Casimir ratio** — symmetry-protected, no dynamics,
exactly positive. The mass ratio `ρ=√(8/3)≈1.633` sits remarkably close to the corpus's physical
excited-glueball threshold `≈1.5 m_g` (a model observation, not a prediction). This is the
sharpest existence proof that **isolation is a real, computable spectral datum** in a confining
theory, and that it degenerates only in the free/gapless limit.

**(W3) 1D chain (quantum rotor) — exact discrete tower `[T]+[V]`.** `E_n=n²/2I`: levels
`0, 0.5, 2, 4.5, 8`, gaps `0.5,1.5,2.5`; C-even tower `0,2,8` ratio `4`. As `I→∞` (free) every fixed
`E_n→0`, giving a continuum from `0` and **no gap** — the gapless limit correctly loses isolation.

**(One-plaquette / single-link `[T-standard]`).** No time direction ⇒ no transfer spectrum; the
character sum is a single time slice whose "tower" is the same `{C₂}` set as (W2) (Laplace
spectrum of the gauge group), so isolation there is again the Casimir ratio, degenerate with 2D
leading order. It fixes the *group-theoretic* content (the ratio is symmetry-fixed) but carries no
dynamical mass; reported as a sanity anchor only.

**Model verdict.** In every solvable *confining* model isolation exists, is exact, and is
symmetry-fixed (Casimir ratio) — it is **not** a small-constant accident. In the free/gapless
limit it correctly vanishes. The hard content is entirely the 4D continuum, not the concept.

---

## 5 · Single-scale covariance: why the *fixed-`δ`* isolation is the gap, and the *ratio* is the datum

The model tests expose a tension. In 2D YM isolation is *free* (Casimir ratio) — yet in 4D it is
the wall. The resolution is **Single-Scale Covariance (SSC)** (witness W2, `G2D §2`):

> **SSC `[T]`.** In a one-scale theory every mass is `c·Λ`, so `Λ★=ρ²m₁²` with `ρ` β-independent,
> and `Λ★−m₁²=(ρ²−1)m₁²`. Hence the *fixed-absolute-`δ`* statement
> `Λ★−m₁²≥δ` **⟺** `m₁²≥δ/(ρ²−1)` **⟺** the gap (given `ρ` bounded).

The witness makes it concrete (`ρ²=8/3`): `m₁²=1→`isolation `1.67`; `m₁²=0.01→`isolation `0.017`.
A **fixed absolute** margin needs `m₁²>0`; a **fixed fractional** margin needs only `ρ>1+c`. Two
honest readings:

- **Clay needs the gap**, i.e. an absolute `δ` in physical units; SSC says that *is* the gap — no
  easier statement hides here. (This is why G2D/STATE-DENSITY both say uniform-`I2` ⟺ gap.)
- **The genuinely new, scale-invariant datum is `ρ≥1+c`** — a dimensionless ratio that is
  acid-clean (the free field has no discrete tower to order, `ρ` undefined/`1`) and, if proven,
  would give uniform `I2` from the *already-proven* strong-coupling gap (SZZ) without touching
  `a→0`. Proving `ρ>1` rigorously is **not** a moment inequality (Thin-gap `[NG]`); it is the
  **order-2** statement of §3.3, and in the strong-coupling polymer it is the same `O(u)`
  cost-vs-rate residue that obstructs the atom (`SUB2-I2-STRONG`: `Cu` vs `Be^{−δt}`).

So: **isolation is "the same structure as the atom" at the *dynamical estimate* level (both are the
one `O(u)` wall, by SSC), but "genuinely independent" at the *spectral-functional* level (different
functional — second support point vs. atom residue; different Karlin tier — order-2 vs order-1).**
Both statements are true and not in tension; the wave states them together.

---

## 6 · The exact missing estimate and its falsification test

### 6.1 Exact missing estimate (E-I2) `[O]`

> **`(E-I2)`.** There exist `c>0` (and `a₀,β₁,L₀`) such that, uniformly for `β` in the physical
> window through `β≈2.5`, `a<a₀`, `L>L₀`,
> ```
>    ρ(L,a,β) := m_{0⁺⁺*}(L,a,β) / m_{0⁺⁺}(L,a,β)  ≥  1 + c ,
> ```
> where `0⁺⁺*` is the **first excited** `0⁺⁺` — the **operative threshold `Λ★`** (not `2m_g`:
> `C/P` cannot separate `gg` from `(gg)(gg)`, so the excited glueball at `≈1.5 m_g` fills below the
> two-body cut, `PROOF-ISOLATION §4.4`). Equivalent forms:
> - **spectral:** `λ₁−λ₂ ≥ (1−e^{−a c m₁}) e^{−a m₁}` on the `0⁺⁺`-projected transfer matrix;
> - **transfer/TP:** an **order-2 total-positivity** (second invariant cone / Karlin flag) of the
>   coarse 0⁺⁺ transfer kernel over `𝒜/𝒢`;
> - **strong-coupling:** the two-exponential remainder `|e^{m₁t}C(t)−f₁²|≤Be^{−δt}`, `δ≈(ρ−1)m₁`,
>   i.e. the `O(u)` self-energy of `λ₁` *resolved against* `λ₂` (`SUB2`) — the same `O(u)` sign as
>   the atom's residue, one tier up.

**Given `(E-I2)`** + the proven SZZ strong-coupling gap `m_{0⁺⁺}≥μ` (and bounded `ρ`) ⇒ **uniform
`I2`** with `δ=(ρ²−1)μ²`; **given** additionally `Λ★>R̄` (weight dominance) ⇒ uniform **atom floor
`f₁²≥r₀`** by Lemma D; **given** the ceiling `χ≤1/λ` ⇒ **mass gap** `m₁²≥r₀λ` (`ch02 eq:sufficient`).

The two side-conditions (`ρ` bounded above, `Λ★>R̄`) are *not* free but are benign in a gapped,
below-threshold-dominated spectrum; they must be carried, not dropped (charter §2).

### 6.2 Falsification test (any valid proof of `I2` must survive all)

- **F1′ free massless `tr F²` — must yield NO atom and `λ₂=λ₁` `[T]+[V]`.** The free massless
  composite has a 2-gluon cut from `0`, `f₁²=0`, `ρ=1` (`W1`). An "isolation theorem" that returns
  `ρ>1` here is **refuted** (it must return `0`, like `H6 §6A` does for the free field).
- **F2′ thin-gap moment-match — must NOT be certified from moments `[T]+[V]` (`W4`).** For every
  finite `K` there is a measure matching `K` moments with isolation `→0`. Any positivity/moment-only
  argument claiming `I2` is **refuted** (decisive, the `I2`-analogue of the thin atom).
- **F3′ double well (Kramers) — no uniform second gap from typical-point `G″>0` `[T]+[V]`.** On
  `V_eff=(g/2)(x²−1)²`, `G″(well)` grows while `λ₁∼e^{−ΔV}→0`; the *second* gap collapses faster.
  A claim of a typical-point order-1 constant controlling `λ₁−λ₂` is **refuted** (`H6 §6B`).
- **F4′ continuum-tower illusion — must not extract `I2` from moments with unknown `K` `[T]+[V]`
  (`W5`).** An argument that "computes `ρ` from moments" must fail when the atom count is unbounded;
  the 8-atom-probe (W5) shows the phantom isolation. Refutes any moment-certificate in 4D.
- **F5′ atomless-gap `:φ²:` — must not confuse "gapped" with "isolated atom" `[T]+[V]`.** Free
  massive `:φ²:` has `inf supp=4m²>0` (gapped) yet `f₁²=0` and `Λ★=m₁²` (no isolation, `W1`). A
  proof that infers `I2` from "gapped" alone is refuted; conversely a valid proof must produce the
  atom *from* the isolation (Lemma A/D), not separately.

---

## 7 · Status ledger

| item | statement | label |
|---|---|---|
| Lemma A | `I2 ⟹ I1` (support gap forces atom) | **[T]** |
| Lemma B | `I1 ⇏ I2` (continuum may touch `m₁²`) | **[T]** |
| Lemma C | `Λ★ ≤ (M₁−f₁²m₁²)/(M₀−f₁²)` (moments upper-bound isolation; dual of Christoffel) | **[T]** |
| Lemma D | `Λ★>R̄ ⇒ f₁²≥M₀(Λ★−R̄)/(Λ★−m₁²)` (isolation ⟹ atom floor) | **[T]** (+side cond.) |
| target `(T★)` = `I2` alone (qualitatively) | atom is downstream, not independent | **[T]** (correction to Q6 §5) |
| thin-gap: moments/positivity cannot lower-bound `I2` | lower end `0`, realised | **[NG]** + **[V]** W4 |
| finite-tower isolation | Prony: `Λ★−m₁²` algebraic in `2K` moments iff `K<∞` | **[T]** + **[V]** W5 |
| order tier | atom = order-1 (Perron); isolation = order-2 (TP2/second cone) | **[T] + [NG]** |
| second cone on `𝒜/𝒢` | absent (orderless quotient) — the wall | **[NG]** (`INVENT-SECOND-CONE`) |
| 2D YM isolation | Casimir ratio `ρ²=8/3`, exact, symmetry-protected | **[T]+[V]** W2 |
| 1D rotor isolation | exact gaps; free limit gapless | **[T]+[V]** W3 |
| free massless `tr F²` | no atom, `ρ=1` — acid PASS | **[T]+[V]** W1 |
| SSC | fixed-`δ` `I2` ⟺ gap; scale-invariant datum = `ρ≥1+c` | **[T]** (+`[C]` on `ρ` bounded) |
| uniform `I2` in 4D `β≈2.5` | `(E-I2)` | **[O]** |
| atom `f₁²≥r₀` in 4D | implied by `I2`+weight-dominance; else P3 at `β<1/48` | **[O]** / **[T]**@`β<1/48` |
| Clay gap | unbuilt | **[O]** · **0%** |

---

## 8 · Sharpest line

**Isolation `I2` is not a corollary of the atom — it is the atom's parent.** A support gap forces the
atom (`I2⟹I1`, Lemma A); the atom does not force a gap (`I1⇏I2`, Lemma B: a continuum may touch
`m₁²` from above); and moments bound the isolation only from **above** (Lemma C, the dual of the
Christoffel bound on `f₁²`), with the lower end **realised at `0`** by the thin-gap measures (W4). So
isolation is at least as hard as the atom, and is a **different, second-level spectral datum**. It
sits one **Karlin tier** higher: the atom is order-1 positivity (Perron/positive kernel), isolation
is order-2 (`λ₁−λ₂`, a second invariant cone) — and that second cone is exactly what the orderless
quotient `𝒜/𝒢` lacks, which is why every order-1 tool (PF, RP, Bushell, BE, Cheeger) is blind to it
and fails as a *tier*, not a constant. In every solvable confining model isolation is **real, exact,
and symmetry-fixed** (2D YM Casimir ratio `ρ²=8/3`; 1D rotor exact gaps), and in the free massless
`tr F²` it is **correctly absent** (acid PASS). In 4D it is `[O]`; by single-scale covariance the
fixed-`δ` form **is** the gap, and the scale-invariant content is the dimensionless ratio
`ρ=m_{0⁺⁺*}/m₁≥1+c` — the same `O(u)` cost-vs-rate residue that walls the atom, one order-2 tier up.
**The missing estimate is `(E-I2): ρ≥1+c` uniformly (order-2 total positivity on `𝒜/𝒢`); the atom
follows from it, not the other way. No QED. Clay 0%.**

---

*Refs.* **Isolation banked:** `unconditional/PROOF-NO-DISSOLUTION-KINEMATIC.md` (`m_{0⁺⁺}<2m_g`,
`[C]`); `PROOF-ISOLATION.md` (operative `Λ`=excited `0⁺⁺*≈1.5m_g`, P3, L1a/L1b);
`PROOF-EDGE-E2.md` (bottom-functional blindness, Temple regress, two-mass estimate);
`G2D-ISOLATION-I2.md` (SSC; vNW pointwise; RP vacuous; `ρ`); `G2F-ATOM-PLUS-ISOLATION.md`
(determinant zero on orbit radius); `STATE-DENSITY-I2.md` (window count `O(1)`, Weyl benign,
`ρ` residual); `SUB2-I2-STRONG.md` (`Cu` vs `Be^{−δt}`; grading grades amplitudes not eigenvalues);
`G2G-TRANSFER-SPECTRAL.md` (Bushell lower bound = `G″` = H1; isolation = `λ₁−λ₂`); `G2H-BLIND-CLEAN.md`
(`(MS)`, moment trap); `SUPPORT-NOT-ATOM.md` (`:φ²:` gapped atomless; support ≠ isolation);
`INEQ-T1-hankel-stieltjes.md` (Christoffel upper bound; atom lower end `0`); `INEQ-T5-sumrule-hunt.md`
(positivity forces minimal residue to zero; Caron-Huot–Pokraka–Zahraee JHEP 01 (2024) 195).
**Logical reduction:** `docs/paper-v2/ch02-squeeze.tex` (squeeze, thin atom, `eq:sufficient`,
`§schur`); `ch05-theorems.tex` (`sec:reslevel` residue–level–gap lemma, `eq:reslevel`);
`ch03-existence.tex` (`§dim`). **Order-2 / second cone:** `INVENT-SECOND-CONE.md`;
`research/v3/H6-HEAD-SPECTRAL.md` §8 (no second cone on `𝒜/𝒢`); `PROOF-GAUSS-INVARIANT-PF.md`
(edge obstruction; free Maxwell); `PROOF-F1SQ.md` (`:φ²:` `f₁²=0`, usc). **Spectral models:**
Migdal 1975 (2D YM solvability); Osterwalder–Seiler, Ann. Phys. **110** (1978) 440 (strong-coupling
gap). **Standard tools as such:** Prony/Padé (finite-tower moment recovery); Kato 1976 VII §3;
von Neumann–Wigner 1929; Källén 1952 / Lehmann 1954; Shen–Zhu–Zhu, CMP **400** (2023) 805,
arXiv:2204.12737 (SZZ gap). **Witness:** `research/v3/w7/r2_isolation.py` (W1–W6, `py -3.12`,
deterministic, numpy/scipy, all PASS). **No third-party downloads; no witness re-run beyond W1–W6;
no `[C]` upgraded to `[T]`; no QED; Clay 0%.**
