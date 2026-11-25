# U6-BINDING-LOGIC — the **exact logical strength** of the binding theorem: a uniform residue floor `w₁≥w₀` implies **only the atom `A`**, is **incomparable with isolation `I`** (both directions witnessed), **does not force the window `ρ²≤C²`** (double well), and proves the **gap `G` only together with the independent ceiling `H1`** (squeeze) — it never proves `I` and never proves `G` alone

> **2026-11-24 · Wave U6 — logical-strength audit of the binding theorem.** The T wave closed on a
> shape for the missing object: *a nonperturbative binding theorem, a uniform lower bound
> `w₁≥w₀>0` on the atom fraction of the minimal four-wall tube, equivalently a leakage bound
> `Σ_{n≥2}w_n≤1−w₀`, strong enough to force the TP₂ order and the window `ρ²≤C²`*
> (`T-WAVE-SYNTHESIS §5`). This note does **not** attack whether such a theorem is true. It fixes what
> such a theorem, *if proved*, would logically buy: does it reach the atom `A`, the isolation `I`, or the
> gap `G`?
> **Reads.** `T-WAVE-SYNTHESIS.md`, `T6-BINDING-ADVERSARIAL.md`, `T3-TP2-ISOLATION.md`,
> `R2-ISOLATION.md`, `R6-ATOM-ADVERSARIAL.md`, `Q6-H1-STATUS.md`; `docs/paper-v2/ch02-squeeze.tex`
> (`eq:squeeze`, `eq:sufficient`, `eq:hstar`, `§schur`), `ch07-nogos.tex` (nested-cone, B2 single phase).
> **Witness:** `research/v3/w7/u6_binding_logic.py` (`py -3.12 -I`, numpy 2.4.6 + scipy 1.18.1,
> deterministic; **27/27 PASS**, exit 0). The double-well sweep reproduces `T6 §6·B` / `R6 §6·C` to
> printed digits.
> **Labels (charter §2).** **[T]** theorem/identity · **[C]** conditional · **[NG]** no-go for a named
> class · **[O]** exact missing step · **[V]** numerical witness. **No QED. Clay 0%.**

---

## 0 · Verdict (up front; no QED)

**The task: the exact logical strength of the binding theorem.** Five statements, each labelled.

1. **Binding implies only the atom, `B ⟹ A`, and nothing more `[T]`.** The residue floor `B` (the
   bottom atom carries a fraction `w₁=f₁²/M₀≥w₀`) is a *weight* datum; it forces `f₁²>0` (that is `A`,
   quantitatively strengthened) and **nothing else**. It does **not** imply isolation, and it does
   **not** imply the gap. This is the headline.

2. **Binding and isolation are logically INCOMPARABLE `[T]` — the T wave's "`w₁≥w₀ ⟺$ isolation`" is a
   conflation and is false as stated.** `B` is a statement about the **weight** of the bottom atom; `I`
   is a statement about the **support/position** of the next atom. Neither implies the other:
   - `B ⇏ I`: a threshold atom (`ν=w₀δ_{m₁²}+` continuum touching `m₁²` from above) has weight but
     `Λ★=m₁²` (`I` false). **[T]+[V]**
   - `I ⇏ B`: a thin isolated atom (`ν=εδ_{μ₁}+(1−ε)δ_{μ₂}`, `μ₂>μ₁`) is isolated (`I` true) but
     `w₁=ε→0` (`B` false), *even with a bounded ratio `ρ²=μ₂/μ₁=2` and the gap `G` holding*.
     **[T]+[V]**
   So `B` and `I` are **two different axes over `A`** (weight vs separation), not one object. The
   T wave's "`w₁≥w₀ ⟺$ leakage `⟺$ isolation `=$ H*`" is correct in *none* of its three equated slots.

3. **Binding does NOT force the window `ρ²≤C²` `[T]+[V]` — the double-well no-go.** The precise form
   "strong enough to force the TP₂ order and `ρ²≤C²`" is **impossible for a weight-only theorem**: the
   Savvidy double well keeps `w₁=O(1)` (witness: `w₁∈[0.962,0.981]` across the barrier sweep `[V]`),
   keeps a large isolation `Δ=μ₂−μ₁` (`I` holds), while `ρ²=μ₂/μ₁→∞` and the gap **Kramers-collapses**
   `λ₁:0.79→4.6×10⁻³` and **`χ→∞` (`H1` fails)**. A residue floor is blind to the level *spacing*; the
   upper window is an independent level-density (`IN-11`) datum, never a consequence of binding. **[T]+[V]**

4. **A proof of the binding theorem does NOT prove the gap; it remains strictly weaker `[T]`, and the
   gap follows only with the independent ceiling `H1`.** Under the standard normalization `M₀=1`
   (`w₁=f₁²`), the squeeze (`ch02 eq:squeeze`, `f₁²/χ≤m₁²≤M₀/χ`) plus a uniform susceptibility ceiling
   `H1: χ≤1/λ` gives
   ```
        B (f₁²≥w₀)  ∧  H1 (χ≤1/λ)   ⟹   m₁² ≥ f₁²/χ ≥ w₀λ > 0   =   G .              [T]
   ```
   This is `ch02 eq:sufficient` exactly. **`H1` is the weakest hypothesis under which binding gives the
   gap**, and it is genuinely independent (a `gap` does not give `B`: the gapped atomless `:φ²:` has
   `G` but `w₁=0`). So the honest answer: **binding proves the atom (`A`), and proves the gap only when
   paired with the orthogonal ceiling `H1`; it never proves isolation `I`.**

5. **Under the bounded-`ρ` window, binding is still NOT equivalent to isolation `[T]`.** The window
   theorem (`T3-B`: `Δ=(ρ²−1)μ₁`, `1+c≤ρ²≤C²`, `Δ≥δ ⟹ μ₁≥δ/(C²−1)`) makes **`I` equivalent to `G`**
   (given the margin), and its "lower bound `ρ²≥1+c`" is **`I` itself** (the second eigenvalue separated
   from the first), **not** the residue floor `B`. The thin isolated atom (bounded `ρ²=2`, `I`, `G`, but
   `w₁→0`) is the witness that the window does **not** deliver `B`. So under the window, `B ⟂ I`; the
   window route to `G` uses `I`, never `B`, and binding adds nothing to it.

**One line.** *A uniform residue floor `w₁≥w₀` is not isolation and not the gap: it implies only the
atom `A`; it is incomparable with the isolation `I` (threshold atom vs thin isolated atom); it cannot
force `ρ²≤C²` (double well: `w₁=O(1)` while `ρ²→∞`, gap→0); and it yields the gap only through the
squeeze together with the independent ceiling `H1`. The T wave's identification of the residue floor
with "isolation = H*" is a weight-vs-support conflation. The binding theorem is a genuine, strictly
weaker sub-target: it proves the atom's residue, never the gap alone, never isolation. No QED. Clay 0%.*

---

## 1 · Definitions and the normalization that matters

`ν≥0` is the Källén–Lehmann measure of the `0⁺⁺` composite `O=tr F²` on `s=m²∈[0,∞)`,
`C(t)=∫e^{−√s\,t}dν(s)`, `f₁²=ν({m₁²})`, `m₁²=inf supp ν`. Moments `M₀=∫dν` (zeroth, total weight),
`M_{−1}=χ=∫s⁻¹dν` (susceptibility). Throughout, **normalize `M₀=1`** (spectral-weight fractions), as
the T wave does; then the absolute residue and the fraction coincide.

| symbol | definition | meaning | axis |
|---|---|---|---|
| **`A`** | `f₁²>0` | atom visible at the bottom of support | weight (positivity) |
| **`B`** | `w₁:=f₁²/M₀≥w₀>0` (uniform in `(L,a)`) | **binding = residue FRACTION floor** | **weight (quantitative)** |
| **`I`** | `Λ★:=inf supp(ν−f₁²δ_{m₁²})>m₁²` | isolation (support gap above the atom) | **position/support** |
| **`G`** | `m₁²>0` | gap (floor of support positive) | position (location) |
| **`H1`** | `χ=∫s⁻¹dν≤1/λ` (uniform) | susceptibility ceiling | moment (`−1`st) |

The T wave's equivalent forms of `B`: **(i) fraction** `f₁²/M₀≥w₀`; **(ii) leakage** `Σ_{n≥2}w_n≤1−w₀`.
*Precision `[T]`.* (ii) equals (i) **only** if (a) the spectrum is purely atomic (no a.c. part) and (b)
`M₀=1`. With a continuum present, `1−w₁=Σ_{n≥2}w_n+M_ac`, so a bound on the *higher-atom* sum alone
(`Σ_{n≥2}w_n≤1−w₀`) does **not** bound `w₁` below — the continuum can carry the weight (witness §6·D:
`w₁=0.005` while `Σ_{n≥2}w_n/M₀=0≤0.9`). The honest leakage is `1−w₁` (include `M_ac`). **A binding
theorem stated as "higher-atom leakage" is strictly weaker than one stated as the bottom fraction.**

*Scheme point `[T]`* (R6-c). Under transport `f₁²↦Z_O²f₁²` (`Z_O²>0` a square), the *fraction*
`f₁²/M₀` and the *ratio* `f₁²/χ` are `Z_O`-invariant (both numerator and denominator scale by `Z_O²`),
while the *absolute* floor `f₁²≥r₀` is not. So `B` (fraction) is the scheme-clean form — but see §4:
the **squeeze-feedable** ratio is `f₁²/χ`, not `f₁²/M₀`.

---

## 2 · The logical lattice: `B ⟹ A`, `B ⟂ I`, `B ⟂ G` (both separations strict)

> **Theorem U6-1 (`B ⟹ A`, strict) `[T]`.**
> `w₁=f₁²/M₀≥w₀>0` with `M₀>0` gives `f₁²≥w₀M₀>0`, i.e. `A`. Strict: `A` holds with `w₁` arbitrarily
> small (no uniform floor — the thin atom, U6-3). So `B` is the *quantitative* strengthening of `A`.
> ∎

> **Theorem U6-2 (`B ⇏ I`; the threshold atom) `[T]+[V]`.**
> `ν=w₀δ_{m₁²}+1_{(m₁²,∞)}ρ` (or a discrete tower `Σ_n w_nδ_{m₁²+1/n}`). Then `w₁=w₀/(w₀+…)` can be
> held `≥w₀` while `Λ★=m₁²` (the a.c./tower touches the atom from above). So `B` holds and `I` fails.
> **Witness §6·A**: `w₁=0.600, A=T, I=F, G=T` (continuum touching); and the finite tower shows the
> isolation **margin has no uniform lower bound** (`1/N: 1.0×10⁻⁴→1.0×10⁻⁵`). ∎

> **Theorem U6-3 (`I ⇏ B`; the thin isolated atom) `[T]+[V]`.**
> `ν=εδ_{μ₁}+(1−ε)δ_{μ₂}`, `μ₂>μ₁`. Then `Λ★=μ₂>μ₁` (`I` holds), `m₁²=μ₁>0` (`G` holds), the ratio
> `ρ²=μ₂/μ₁` is **bounded** (`=2` for `μ=(1,2)`), yet `w₁=ε→0` (`B` fails). So `I` — even with a
> **bounded window and the gap** — does not give `B`. **Witness §6·B**: `ε=10⁻⁶, I=T, G=T, ρ²=2,
> w₁=10⁻⁶`. ∎

> **Corollary (U6-2 ∧ U6-3): `B` and `I` are INCOMPARABLE `[T]`.** Neither implies the other; both
> separately imply `A` (`B⟹A`; `I⟹B_spec⟹A`, T6-1). They are **two different strengthenings of `A`**:
> `B` on the *weight* axis, `I` on the *support* axis. **This is the precise correction of the T wave's
> "`w₁≥w₀ ⟺$ isolation`": they are incomparable coordinates of `ν`, not one object.** ∎

> **Theorem U6-4 (`B ⇏ G`; massless atom and double well) `[T]+[V]`.**
> *Measure level.* `ν=w₀δ_0+(1−w₀)δ_M`: `w₁=w₀` (`B`), `Λ★=M>0` (`I`), but `m₁²=0` (`G` fails) and
> `χ=∞` (`H1` fails). Witness §6·C: `w₁=0.5, I=T, G=F, χ=∞`. *Physical level.* The Savvidy double
> well (§3) has `B` and `I` with `G` collapsing. So `B ⇏ G`. ∎

> **Theorem U6-5 (`G ⇏ B`; the gapped atomless) `[T]+[V]`.**
> `ν=1_{[4m²,∞)}ρ` (free massive `:φ²:`): `m₁²=4m²>0` (`G`), but `f₁²=0` (`A`, `B` fail). So `G` and
> `B` are independent. Witness §6·B: `m₁²=4, G=T, w₁=0`. ∎

**The lattice.** Three axes over the atom: **weight** (`A`→`B`), **support** (`I`), **location** (`G`).
```
        B ─(strict)─► A ◄─(strict)─ I            (B = quantitative A ; I = separation of A)
        B  ⟂  I          (neither implies the other; U6-2, U6-3)
        G  ⟂  B          (gapped atomless ⇏ B ; massless/double-well B ⇏ G ; U6-4, U6-5)
        G  ⟂  A         (gapped atomless ; I ⟹ A but G alone does not)
```
The T wave's "true axes `A, I`, gap orthogonal" is right about `A, I, G`; the error is calling `B`
"isolation." `B` is **not** a fourth axis — it is the quantitative version of `A`, sitting on the
*weight* line, **orthogonal to `I`**.

---

## 3 · Binding does not force the window: the double-well no-go `[T]+[V]`

The T wave wants the binding theorem "strong enough to force the TP₂ order and the window `ρ²≤C²`."
The double well proves a **weight-only** theorem cannot.

Witness §6·E (SUSY/transfer operator `H=−∂²+W`, `W=V'²/4−V''/2`, `V=(g/2)(x²−1)²`, carrier `O=x`):

```
   g      λ₁=GAP       w₁=f₁²/M₀    χ        ρ²=λ₂/λ₁     Δ=λ₂−λ₁
  0.5   7.92e-01      0.9737      1.29       8.6          6.05
  2.0   7.49e-01      0.9619      1.07      15.6          10.9
  8.0   1.19e-01      0.9636      7.43     185           21.9
 12.0   2.48e-02      0.9740     37.2     1257           30.8
 16.0   4.58e-03      0.9812    207       9481           43.4
```

The **fraction `w₁` stays `≈0.97`** (`B` holds, uniformly across the sweep), the **isolation `Δ` stays
large** (`I` holds), while the **gap `λ₁` Kramers-collapses** (`G` fails), **`χ→∞`** (`H1` fails), and
**`ρ²→∞`** (the upper window fails). Hence:

> **No-go U6-NG-W (`B ⇏ ρ²≤C²`) `[T]/[NG]`.** In the class of weight-data (functionals of the atom
> weights `{w_n}`), **no residue floor bounds the ratio `ρ²=μ₂/μ₁` from above, nor bounds `χ`, nor
> bounds `m₁²` from below.** The double well has `w₁=O(1)` and `ρ²→∞`, `χ→∞`, `m₁²→0`. The upper
> window `ρ²≤C²` is an independent level-density (`IN-11`) datum; the susceptibility ceiling `H1` is an
> independent moment datum. **"Binding forces the window" is false; binding is blind to both.** ∎

```
       weight floor  w₁≥w₀   ──╳──►   spacing  ρ²≤C²        (U6-NG-W: double well)
                              ──╳──►   ceiling  χ≤1/λ        (double well: χ→∞)
                              ──╳──►   isolation I          (U6-2: threshold atom)
```

---

## 4 · The two routes to the gap, and where binding actually enters `[T]`

There are **two independent sufficient routes** to `G`, and they consume *different* data:

**(Route 1 — window / position; `T3-B`, no binding).** `I` (`ρ²≥1+c`) + `1+c≤ρ²≤C²` + margin
`Δ=(ρ²−1)μ₁≥δ ⟹ μ₁≥δ/(C²−1)>0`. Pure position/eigenvalue data; **`B` is absent and unnecessary.**
In bounded-`ρ` models (rotor `ρ²=4`, Mehler `2`, box `4`, 2D YM `8/3` — witness §6·F) this gives a
positive gap floor.

**(Route 2 — squeeze / weight; `ch02 eq:sufficient`, needs the ceiling).** `B` (residue floor) +
`H1` (`χ≤1/λ`) `⟹ m₁²≥f₁²/χ≥w₀λ>0`. Witness §6·G: `f₁²=0.5, χ=0.5417≤1/1.6, m₁²=1≥f₁²/χ=0.923≥w₀λ=0.8`.
**Binding enters ONLY here, and only with `H1`.**

> **Theorem U6-6 (the gap-feeding datum) `[T]`.**
> The squeeze feeds on the **ratio `f₁²/χ`** (residue-to-susceptibility), not on the **fraction
> `f₁²/M₀`**. Precisely: `f₁²/χ≥c ⟹ m₁²≥c` already (squeeze lower rung, no extra input); the
> `Z_O`-invariant ratio `f₁²/χ` is the cleanest sufficient datum. The T wave's `w₁=f₁²/M₀` is a
> *different* ratio (`M₀` vs `χ`); with `M₀=1` fixed and `H1` supplying `χ≤1/λ`, one has
> `f₁²/χ=w₁/χ≥w₀λ`, recovering Route 2. Without either the `M₀`-normalization or `H1`, the fraction
> form does **not** alimenter the squeeze. ∎

**Consequence for "does binding prove the gap?":**
- **Under the bounded-`ρ` window** (`I` holds, `ρ` bounded, margin): the window already gives `G`;
  **`B` contributes nothing** (Route 1). And `B` is not even implied (§2, thin isolated atom). So in
  the regime the T wave wants, binding is **superfluous to the gap**.
- **Outside the window** (`ρ` unbounded, or `I` fails): `B` alone cannot give `G` (double well,
  massless atom). `B+H1` does (Route 2) — but that is the squeeze, and `H1` is the real ceiling.
- **Never** does `B` give `I`, and never `G` alone.

So the **exact strength**: *`B` is sufficient for `A`; `B+H1` is sufficient for `G`; `B` is necessary
for neither `I` nor `G`; and `B` is not needed for the window.*

---

## 5 · The three separation models (the task's explicit battery)

| separation the task asks | witness | data | label |
|---|---|---|---|
| **`w₁>0` but not `G`** | **Savvidy double well** (and massless atom) | `w₁≈0.97` uniform, `λ₁→4.6×10⁻³→0`, `ρ²→∞`, `χ→∞`; massless: `w₁=0.5, m₁=0, χ=∞` | **[T]+[V]** |
| **`I` but not `G`** | **bounded vs unbounded `ρ`** | bounded (`ρ²=2,4,8/3`) ⇒ `G` (rotor/Mehler/box/2D); unbounded (`ρ²→∞`) ⇒ `I` with `G→0` (double well; also the massless isolated atom `I ∧ ¬G`) | **[T]+[V]** |
| **`G` but not `w₁`** | **gapped atomless `:φ²:`** | `m₁²=4m²>0` (`G`), `f₁²=0` (`w₁=0`) | **[T]+[V]** |
| *(extra)* **`B` but not `I`** | threshold atom / tower | `w₁=0.6`, `Λ★=m₁²` | **[T]+[V]** |
| *(extra)* **`I` but not `B`** | thin isolated atom | `ε=10⁻⁶`, `I`, `G`, `ρ²=2`, `w₁→0` | **[T]+[V]** |

The **double well does double duty**: it is *simultaneously* the "`w₁>0` but not `G`" witness **and**
the "`I` but not `G` under unbounded `ρ`" witness — exactly because `B` and `I` are both robust while
the *position* (gap) and the *ratio* (window) are the fine-tuned objects. This is the structural
statement of the whole corpus (residue/isolation robust, gap fine-tuned), now sharpened: **the gap is
the only coordinate that the double well kills, and neither `B` nor `I` survives as a proxy for it
once `ρ` is unbounded.**

---

## 6 · Witness — `research/v3/w7/u6_binding_logic.py` (`py -3.12 -I`, 27/27 PASS)

Deterministic, numpy 2.4.6 + scipy 1.18.1, no MC, ~3 s. Quoted output (abridged):

```
PART 1  B vs I vs G at the measure level:  B is INCOMPARABLE with I
(1a)  tower w1=0.848 A=T I=T ; iso-margin 1.0e-04 -> 1.0e-05  (NO uniform isolation)
(1a') threshold-touch: w1=0.600 A=T I=F G=T   => B & ~I exists (B =/=> I)
(1b)  thin eps=1e-6: w1=1e-6  I=T G=T chi=0.500 rho2=2.00  => I & ~B, bounded window, G
(1c)  gapped atomless: m1^2=4 G=T w1=0   => G & ~B (G =/=> B)
(1d)  massless: w1=0.500 I=T G=F chi=inf  => B & ~G (B =/=> G)

PART 2  double well: B & I robust, GAP collapses, rho^2 -> inf
   g     lam1=GAP      w1        chi       rho^2
  0.5  7.921e-01   0.9737    1.29e+00      8.6
 16.0  4.575e-03   0.9812    2.07e+02   9481.4
  => w1 O(1) & I large while GAP collapses, chi->inf, rho^2->inf   (B & I not enough for G)

PART 3  window:  I + (1+c<=rho^2<=C^2) + margin => G ;  B absent
  rotor rho^2=4, Mehler 2, box 4, 2D SU(2) 8/3  => positive gap floor
  thin isolated atom: I + rho^2=2 + G  but w1=1e-6  => window does NOT give B
  unbounded rho^2 at fixed Delta => mu1=Delta/(rho^2-1)->0 (double-well signature)

PART 4  squeeze:  B_abs + H1 => G ;  B alone =/=> G
  f1^2=0.50 chi=0.5417 H1(<=0.625)=T => m1^2=1.0 >= f1^2/chi=0.923 >= r0*lam=0.800

PART 5  leakage:  sum_{n>=2}w_n<=1-w0 holds while w1=0.005 (continuum carries weight)
         => higher-atom leakage is NOT equivalent to the bottom fraction

PART 6  battery: B=>A (0 viol), I=>A (0 viol);  B&~I, B&~G, I&~B all constructed.
SUMMARY: 27 checks, 0 FAIL   ALL PASS
```

Witness, not proof; the theorems are §2–§4, the numerics are the structural illustrations. The
double-well sweep reproduces `T6 §6·B` / `R6 §6·C` (same `λ₁` Kramers collapse, same large `Δ`,
same `ρ²→∞`); the `w₁=f₁²/⟨x²⟩` column is new and shows the **fraction** stays `O(1)` while `χ` (a
different moment) diverges — the precise fact that `B` (weight) and `H1` (susceptibility) are distinct.

---

## 7 · Strongest no-go, and the weakest hypothesis

### 7.1 The strongest no-go

> **U6-NG1 (weight-independence no-go) `[T]/[NG]`.** *Class = residue-floor / weight-only functionals
> of the Källén–Lehmann weights `{w_n}`.* In this class the binding predicate `B` is **logically
> independent of both isolation `I` and gap `G`**:
> - `B ⇏ I` (threshold atom, tower: weight without support gap),
> - `I ⇏ B` (thin isolated atom: support gap without weight — even under a bounded window and a gap),
> - `B ⇏ G` (massless atom; double well: weight with `ρ²→∞` and gap→0),
> - `G ⇏ B` (gapped atomless `:φ²:`).
> Therefore **no biconditional** "`B ⟺ I`", "`B ⟺ G`", or "`B ⟺ I`-under-the-window" holds, and the
> T wave's "`w₁≥w₀ ⟺$ leakage `⟺$ isolation `=$ H*`" is a **weight-vs-support conflation, false as
> stated**. The only implications are `B⟹A` and `(B∧H1)⟹G`. A weight-only binding theorem reaches
> **at most `A`**, and the gap **only** with the independent ceiling `H1`. ∎

The sharpest *model-level* statement: **the double well is a `B∧I`-without-`G` object with `ρ²→∞`,
and the thin isolated atom is an `I`-with-`G`-without-`B` object with bounded `ρ²`** — together they
prove `B` neither forces nor follows from the isolation/gap coordinates.

### 7.2 The weakest hypothesis under which binding gives the gap

> **U6-W (weakest gap hypothesis) `[T]`.** Under the standard normalization `M₀=1`, the binding
> theorem gives the gap **iff** the susceptibility is uniformly bounded:
> ```
>    B (f₁²≥w₀)  ∧  H1 (χ≤1/λ uniform)   ⇒   m₁² ≥ w₀λ > 0  (= G).              [T]
> ```
> `H1` is **necessary and sufficient** for the squeeze route: necessity is the massless/thin cases
> (`χ=∞ ⇒` no floor), sufficiency is `ch02 eq:squeeze`. There is **no hypothesis** making `B` alone
> give `G` (`U6-NG1`). Equivalently, the cleanest `Z_O`-invariant sufficient datum is the **ratio**
> `f₁²/χ≥c` (which gives `G` directly, U6-6). **The weakest hypothesis is `H1` (or the equivalent ratio
> floor `f₁²/χ≥c`); binding alone is strictly weaker.** ∎

**Reading of the whole U6 wave.** The binding theorem the T wave wants is a **genuine, strictly
weaker sub-target**: it certifies the atom's *residue* (`A`, quantitatively), and — only when joined to
the independent ceiling `H1` — the gap. It does **not** certify isolation, does **not** force the
window, and is **not** implied by the gap. The campaign's "`Clay = H1 + atom`" pair remains the correct
*skeleton*; U6 fixes the atom's slot: the binding theorem fills the **residue** half of that pair
(`f₁²≥r₀`), and the **ceiling** `H1` is the other, independent half. "Binding ⟹ gap" without `H1` is
false in every model class, and "binding ⟺ isolation" is false in both directions.

---

## 8 · Falsification battery for any proposed binding theorem `[T]/[V]`

A proposed theorem claiming, uniformly in `(L,a)` at `β=2.5`, a residue floor `w₁≥w₀>0` (and perhaps
claiming `I` or `G`) must survive **all eight**; failing any refutes it. Each is an exact /
machine-checked countermodel in the corpus or this note.

- **U6-FB1 — free field / `:φ²:` acid test (must return `f₁²=0`) `[T]+[V]`.** The free massive composite
  has `w₁=0`; the massless free field `f₁²=0`, `χ=∞`. A "binding theorem" that fires here is **vacuous**
  (Schur-concave obstruction, `ch02 §schur`). *Acid test.*
- **U6-FB2 — threshold / BIC (must NOT claim isolation) `[T]+[V]`.** Feed the threshold atom
  (`w₁=0.6`, `Λ★=m₁²`). A *weight* theorem firing here is **correct** (weight exists) — but a theorem
  that claims `I` (or "binding ⇒ isolation") is **refuted** (`B ∧ ¬I`).
- **U6-FB3 — double well (must NOT claim `G` or `ρ²≤C²`) `[T]+[V]`.** Feed `V=(g/2)(x²−1)²`:
  `w₁=O(1)`, `Δ` large, but `λ₁→0`, `χ→∞`, `ρ²→∞`. A theorem claiming the **gap**, a **uniform
  `ρ²≤C²`**, or a **uniform `χ` bound** is refuted. *Refutes "binding forces the window/gap".*
- **U6-FB4 — thin isolated atom (must return `w₁→0`; must NOT infer `B` from `I`) `[T]+[V]`.** Feed
  `ν=εδ_{μ₁}+(1−ε)δ_{μ₂}`, bounded `ρ²`, `G`. A theorem that "proves `B` from isolation / the window"
  is refuted (`I ∧ G ∧ bounded-ρ` with `w₁→0`). Also refutes **uniformity**: a thinning sequence
  `w₁(L,a)→0` satisfies pointwise `w₁>0` but not uniform `≥w₀`.
- **U6-FB5 — massless atom (must NOT claim `G` from `B`) `[T]`.** `ν=w₀δ_0+(1−w₀)δ_M`: `B`, `I`,
  but `m₁²=0`, `χ=∞`. A theorem concluding `G` from `B` alone is refuted.
- **U6-FB6 — leakage precision (must bound `1−w₁`, the bottom fraction) `[T]+[V]`.** A theorem bounding
  only the **higher-atom** sum `Σ_{n≥2}w_n` is refuted by a continuum-dominated measure (`w₁=0.005`
  while `Σ_{n≥2}w_n/M₀=0`). The honest bound is on `1−w₁` (include the a.c. mass) or directly on `f₁²`.
- **U6-FB7 — scheme / `Z_O` / dimension (floor must be a ratio) `[T]`.** Under `f₁²↦Z_O²f₁²` the
  absolute floor `r₀` rescales (`Z_O²`); only the sign and the **ratio** `f₁²/M₀` or `f₁²/χ` are
  invariant. A theorem producing an **absolute** `r₀` without `Z_O²`-covariance or `a²`-scaling is
  scheme/dimension-wrong. *The honest floor is a ratio* (R6-c, Q6-F5).
- **U6-FB8 — carrier + channel specificity `[T]+[V]`.** On the 2D **holonomy** carrier `ρ²=8/3`; on the
  **`tr F²`** carrier the same 2D theory has `f₁²=0` (`tr F²=H` conserved). A **carrier-blind** theorem
  wrongly assigns the holonomy's binding to `tr F²` (T6 §Part4). The theorem must name `O` and the channel.

**What a valid binding theorem therefore is.** A **carrier-and-channel-specific**, `Z_O`-covariant,
`a²`-scaled, **uniform** proof of the bottom residue `w₁=f₁²/M₀≥w₀` that: returns `0` on the free
composite (FB1), fires correctly on the threshold atom without claiming `I` (FB2), does **not** claim
`G` or the window on the double well (FB3), does **not** claim `B` follows from `I` (FB4), does **not**
claim `G` from `B` alone (FB5), bounds the bottom fraction (FB6), is a ratio (FB7), and is operator-
labelled (FB8). **Such a theorem proves `B` (hence `A`), and — with `H1` — the gap; it proves neither
`I` nor `G` alone.** It is exactly what `P3` does at `β<1/48` (a residue floor, no isolation), one
scheme-clean, uniform tier up.

---

## 9 · Status ledger

| item | statement | label |
|---|---|---|
| **U6-1** `B ⟹ A` (strict) | residue fraction floor forces atom visibility | **[T]** |
| **U6-2** `B ⇏ I` | threshold atom / tower: weight without support gap | **[T]+[V]** |
| **U6-3** `I ⇏ B` | thin isolated atom: gap + bounded window without weight | **[T]+[V]** |
| **Corollary** `B ⟂ I` (incomparable) | the T wave's "`B ⟺ I`" is a weight-vs-support conflation | **[T]** (correction) |
| **U6-4** `B ⇏ G` | massless atom; double well `w₁=O(1)`, `G→0` | **[T]+[V]** |
| **U6-5** `G ⇏ B` | gapped atomless `:φ²:` (`G`, `w₁=0`) | **[T]+[V]** |
| **U6-NG-W** `B ⇏ ρ²≤C²`, `B ⇏ χ≤1/λ` | double well: `w₁=O(1)` while `ρ²→∞`, `χ→∞` | **[T]/[NG]** |
| **U6-6** squeeze feeds on `f₁²/χ`, not `f₁²/M₀` | the gap-feeding ratio | **[T]** |
| **U6-NG1** weight-independence no-go | `B` independent of `I,G`; only `B⟹A`, `(B∧H1)⟹G` | **[T]/[NG]** |
| **U6-W** weakest gap hypothesis | `B ∧ H1 ⟹ G`; `H1` necessary+sufficient (or `f₁²/χ≥c`) | **[T]** |
| **leakage precision** | `Σ_{n≥2}w_n≤1−w₀` ⟺ `w₁≥w₀` only pure-atomic & normalized | **[T]+[V]** |
| **window lower bound `ρ²≥1+c`** | IS isolation `I` (eigenvalue ratio), **not** the residue floor `B` | **[T]** (correction) |
| **double well** `B∧I∧¬G`, `ρ²→∞` | the canonical weight-without-gap separator | **[T]+[V]** |
| **thin isolated atom** `I∧G∧bounded-ρ∧¬B` | the canonical separation-without-weight separator | **[T]+[V]** |
| **gapped atomless** `G∧¬B` | free `:φ²:` — gap without residue | **[T]+[V]** |
| **binding theorem, if proved** | proves `A`; gap only with `H1`; **not `I`**, **not `G` alone** | **[T]** (conditional on the theorem) |
| **P3 at `β<1/48`** | an unconditional residue floor `f₁²≥r₀`, **no `I`** — a `B` (cell (a)) anchor | **[T]** |
| **Clay mass gap** | unbuilt | **[O] · 0%** |

---

## 10 · Sharpest line

**The binding theorem is not isolation and not the gap.** A uniform residue floor `w₁=f₁²/M₀≥w₀` is a
statement about the **weight** of the bottom atom; isolation `I` is a statement about the **position** of
the next atom; the gap `G` is a statement about the **location** of the bottom. They are three different
coordinates of the Källén–Lehmann measure, and a weight floor forces none of the position ones:
`B⟹A` (strict) but `B⇏I` (threshold atom: weight with continuum touching) and `I⇏B` (thin isolated
atom: `I` and `G` and bounded `ρ²=2` with `w₁→0`), so `B` and `I` are **incomparable** — the T wave's
"`w₁≥w₀ ⟺$ isolation`" is a weight-vs-support conflation, false in both directions. The double well is
the decisive witness that **no weight-only theorem can force the window or the gap**: it keeps
`w₁≈0.97` and a large `Δ` while `ρ²→∞`, `χ→∞`, and the gap Kramers-collapses to zero. The only way a
residue floor becomes the gap is the **squeeze with the independent ceiling `H1`**
(`B∧H1 ⟹ m₁²≥f₁²/χ≥w₀λ>0`, `ch02 eq:sufficient`); the clean `Z_O`-invariant form is the ratio
`f₁²/χ≥c`, which gives the gap directly. Under the bounded-`ρ` window, the gap comes from **`I`**
(the window's lower bound `ρ²≥1+c` *is* isolation), and `B` is both unnecessary and un-implied. So the
exact strength is: **binding proves the atom `A`; `binding + H1` proves the gap `G`; binding proves
`neither I nor G alone`, is not implied by `G`, and cannot force the window.** It is a genuine,
strictly weaker sub-target — the residue half of `Clay = H1 + atom` — not the gap, not the isolation.
**No QED. Clay 0%.**

---

*Refs.* **Binding shape / targets:** `T-WAVE-SYNTHESIS.md` §5 (binding theorem `w₁≥w₀`, leakage
`Σ_{n≥2}w_n≤1−w₀`); `T2-TUBE-OVERLAP.md` (`w₁=f₁²/M₀`, leakage = "isolation = H*" — corrected here);
`T6-BINDING-ADVERSARIAL.md` (`I⟹B_spec⟹A`, four-cell lattice, `FB1–FB6`); `T3-TP2-ISOLATION.md`
(window `T3-B`, "lower `ρ²≥1+c` = binding half" — corrected here to "= isolation"); `R2-ISOLATION.md`
(Lemma A `I⟹A`, Lemma B `A⇏I`, thin-gap, SSC); `R6-ATOM-ADVERSARIAL.md` (double-well residue robust /
gap collapses, `R6-a` gapped atomless, `R6-c` `Z_O`); `Q6-H1-STATUS.md` (H1 strictly weaker, thin
atom, `F1–F5`). **Paper:** `docs/paper-v2/ch02-squeeze.tex` (`eq:squeeze`, `eq:sufficient`, `eq:hstar`,
`§schur`, thin-atom Prop.); `ch07-nogos.tex` (Schur master, nested-cone, B2 single phase). **Strong
coupling:** `unconditional/PROOF-P3-UNCONDITIONAL.md` (residue floor, no `I₂`). **Witness:**
`research/v3/w7/u6_binding_logic.py` (`py -3.12 -I`, numpy+scipy, deterministic, 27/27 PASS,
double-well reproduces `T6 §6·B` / `R6 §6·C`). **No third-party downloads; no `[C]` upgraded to `[T]`;
no QED; Clay 0%.**
