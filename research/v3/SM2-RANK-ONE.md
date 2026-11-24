# SM2-RANK-ONE — the rank-one soft-mode theorem attempt: exact overlap lemma `[T]`, exact margin identity `[T]`, refined finite-`c` bound `[T]/[V]`, but an exact **ceiling** `m_soft ≤ G″` that forbids an O(1) soft-mode margin `[T]`; head is H1 via the BL term `[O]`

> 2026-11-24 · Attack surface **SM2**: the Savvidy soft mode, the last open mode of the UV face.
> Soft mode = the transverse-magnetic (TM) Savvidy **negative** mode `φ₋`; the defect is **rank one**.
> Task: a positive lower bound on the **soft-mode margin** `m_soft := S − D_Sav` at `β=2.5`, from the
> exact rank-one structure.
> Reads: `H2-DEFECT-BOUND.md` (`D=D₊+D_Sav`, `D_Sav=E|⟨φ₋|c⟩|²/(gB)`), `SAVVIDY-DTERM.md`,
> `SAVVIDY-NEGMODE-EXACT.md` (`Var(φ²)~2/λ₄`), `SAVVIDY-UNIFORM-SOFTENING.md`, `SAVVIDY-2LOOP.md`,
> `H1-RG-DESCENT.md` (head = H1, soft-mode failure).
> Witness: `research/v3/w7/sm2_rank_one.py` (`py -3.12 -I`, numpy only; output quoted §5).
>
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional/model · **[O]** exact
> missing step · **[NG]** no-go for a named class · **[V-num]** this-round numerics.

---

## 0 · Verdict (no QED; Clay 0%)

The rank-one structure is the right coordinate and yields three **exact** results — and one **exact
obstruction** that decides the task:

1. **Overlap lemma `[T]`.** `⟨φ₋|c⟩ = ⟨H_cf φ₋|v⟩` — a *single* coarse–fibre Hessian matrix element
   (verified to `10⁻¹⁵`, §5 (1)). It is **O(1) at the head** `k~√(gB)`, and **suppressed** by three
   named mechanisms off the head: **(i) scale separation** (deep IR, `k≪√(gB)`, the transverse
   coupling `∝k` ⇒ `|overlap|²∝k²→0`); **(ii) one-octave support** (`|p_⊥|~√(gB)`, no suppression
   *in band*, O(1)); **(iii) helicity/gauge projection** (only the `+1` TM projection of `c` couples;
   a longitudinal / pure-gauge `v` has *zero* overlap). So the head is the **worst case**, `O(1)`.

2. **Exact margin identity `[T]`.** `G″ = S − D = (S − D_Sav) + |D₊| = m_soft + |D₊|`, where
   ```
   m_soft := S − D_Sav = E[ h_n − ⟨K₊⁺c,c⟩ ]
   ```
   is the curvature with the negative mode **excised**, and `D₊ ≤ 0` is the Brascamp–Lieb buy-in
   (sign only). The precise inequality for a positive margin `δ` is
   ```
   m_soft ≥ δ  ⟺  E[h_n] − E⟨K₊⁺c,c⟩ − E|⟨φ₋|c⟩|²/(gB) ≥ δ ,     (★)
   ```
   with constants `gB` (negative-mode gap), `λ₄` (quartic), and the measure weight `ρ`.

3. **Refined rank-one bound `[T]/[V-num]`.** `D_Sav ≤ g² Var(gB,λ₄)` with the **exact finite-`c`**
   `Var(φ²)`, *not* the `c→∞` asymptote `2/λ₄`. **Correction to the corpus:** `Var(φ²)≤2/λ₄` is
   **false** as a strict bound — the exact quadrature peaks at `2.21/λ₄` (`gB≈4`, `+10%`). At the
   **head** (`gB=O(1)`, small-`c`) the exact value is only `0.2–0.5 ×(2/λ₄)`, i.e. **2–3× tighter**
   than `2g²/λ₄`. With this, `m_soft ≥ S_bare − g²Var(gB,λ₄)` is **positive** for natural-to-stiff
   `λ₄` and moderate `gB` (positivity region §3).

4. **The decisive exact obstruction `[T]`.** Because `D₊ ≤ 0` has *no lower bound* (BL gives the sign,
   not the magnitude), the margin is **capped above**:
   ```
   m_soft = G″ − |D₊|  ≤  G″ .     (★★)
   ```
   The measured `G″ ≈ +0.049 = 3.5%·E[S]` (`L=8,β=2.4`) therefore **caps any positive soft-mode margin
   at `3.5%·S`**, and `m_soft` can be *arbitrarily negative*. So the observed positivity is
   **BL-carried, not rank-one-carried**: proving `m_soft>0` is *sufficient* for `G″>0` but is
   **strictly stronger than needed and not supported by the data**. **The rank-one soft-mode route
   cannot deliver an O(1) margin; at best a thin `O(3%)·S` one, and possibly none.**

   *(Consistency pins `gB`.)* The bound `m_soft ≥ m_toy=S−g²Var` and the ceiling `m_soft ≤ G″` must
   coexist, i.e. `m_toy ≤ G″` ⇔ `g²Var(gB,λ₄) ≥ S−G″ ≈ 0.965·S`. At natural `λ₄=1.6` this requires
   **`gB ≳ 2.1`** — the marginal regime where `m_toy` is *thin* (`≲0.05`, dropping to `≤0` as `gB`
   grows). The **`gB=1` corner (`m_toy=0.55`) is excluded by the data**: the measured `D_Sav≈S` forces
   the head into the thin regime, and `(★★)` confirms it. So the two exact statements together **pin
   the head's negative-mode curvature to `gB≳2.1`**, exactly where a positive margin (if any) is
   `O(3%)·S`. Witness (7c).

**Bottom line.** The rank-one structure gives the exact overlap, the exact margin, and a *conditional*
positive margin (computable region), but an exact ceiling `m_soft ≤ G″` plus the `O(1)` head overlap
deny a uniform positive `δ`. The head remains **H1**, now shown to be carried by the **BL term** `|D₊|`
whose *magnitude* (not just sign) is the missing input. **Clay 0%; no QED.**

---

## 1 · The exact overlap `⟨φ₋|c⟩` — O(1) at the head, suppressed off it `[T]`

**Setup.** `K = K₊ − gB φ₋φ₋ᵀ`, `φ₋ ∈ ker K₊`, `‖φ₋‖=1` (TM lowest Landau level, `λ₋=−gB`).
Coarse mode `v` (`‖v‖=1`), fibre coupling `c = Π₀ H v = H_fc v` (`H=Hess V`, `Π₀`=fibre projector),
overlap `⟨φ₋|c⟩ = φ₋ᵀH_fc v`.

### 1.1 The overlap lemma `[T]`

Since `φ₋` is a fibre vector (`Π₀φ₋=φ₋`) and `H` symmetric,
```
⟨φ₋|c⟩ = φ₋ᵀH_fc v = φ₋ᵀ Π₀ H v = φ₋ᵀ H v = (H_cf φ₋)ᵀ v = ⟨H_cf φ₋|v⟩ .
```
So the overlap is **exactly the coarse–fibre Hessian matrix element** `⟨H_cf φ₋|v⟩` — the vector `v`
dotted with the image `H_cf φ₋` of the negative mode under the coarse–fibre block. Witness (1):
`⟨φ₋|c⟩ = ⟨H_cf φ₋|v⟩ = 0.98221842`, error `0.0`.

This is the clean algebraic object behind `D_Sav`: **one rank-one matrix element of the Hessian
between the TM negative mode and the coarse mode.**

### 1.2 O(1) vs suppressed — the three mechanisms `[T in the LLL model / structural]`

`H_cf` is (covariantly) translation-covariant; the matrix element is supported when the coarse
momentum `k` can reach the magnetic octave of `φ₋`. In the LLL model the squared overlap factorizes
(witness (2)):
```
|⟨φ₋|c_k⟩|²  ~  [ transversality: k²/(k²+gB) ] × [ LLL UV decay: e^{−k²/(gB)} ] ,   r ≡ k/√(gB).
```
| `k/√(gB)` | `|⟨φ₋|c_k⟩|²` | regime / mechanism |
|---|---|---|
| `0.0–0.1` | `≈0` | **scale separation** — transverse coupling `∝k` ⇒ `∝k²→0` |
| `0.3–1.4` | `0.08–0.20` | **one-octave support** — O(1), *peak 0.21 at `r≈0.79`* |
| `≥2` | `≲0.015` | **LLL UV decay** `e^{−k²/gB}` |

Three named mechanisms, as requested:
- **(i) Scale separation.** Deep IR `k≪√(gB)`: the TM negative mode and a long-wavelength *transverse*
  coarse mode decouple through the derivative/transversality of the cubic coupling, `|overlap|²∝k²→0`.
  *(The Gaussian LLL alone peaks at `k=0`; the decay is the transverse-coupling factor `k²/(k²+gB)`,
  not the Gaussian — stated honestly.)*
- **(ii) One-octave support.** In band `|k|~√(gB)` there is **no suppression** (`k²/(k²+gB) ~ O(1)`,
  `e^{−r²} ~ O(1)`); the overlap is **O(1)**. This is the head, and it is the *worst* case.
- **(iii) Helicity / gauge projection.** `φ₋` is the `+1`-helicity (TM) state. Only the `+1`-helicity
  projection of `c` contributes: `⟨φ₋|c⟩ = ⟨φ₋|P_{(+1)}c⟩`. A **longitudinal / pure-gauge** coarse
  mode `v` (the *softest* `k→0` mode) has `P_{(+1)}c = 0` ⇒ **zero overlap**. So the gauge-projected
  soft mode is invisible to the negative mode; the overlap lives in the *transverse-magnetic* channel.

**Conclusion.** The overlap is **O(1) at the head** (`k~√(gB)`), suppressed only off the head
(deep-IR by scale separation, UV by LLL decay, gauge sector exactly by helicity projection). The head
is where `D_Sav` is largest and where the rank-one bound is loosest.

---

## 2 · The exact margin identity, the inequality, and its constants `[T]`

### 2.1 The margin identity `[T]`

From `D = D₊ + D_Sav`, `D₊ ≤ 0` (BL on the `K₊⪰0` positive fibre), and `G″ = S − D`:
```
G″ = S − D₊ − D_Sav = (S − D_Sav) + (−D₊) = m_soft + |D₊| ,    m_soft := S − D_Sav .
```
Equivalently (pointwise, before averaging)
```
m_soft(x) = h_n(x) − ⟨K₊⁺c,c⟩(x) = (the Gaussian marginal stiffness with the negative mode excised).
```
Witness (1): on a random rank-one cell, `G″=(S−D_Sav)+|D₊|` matches `S−(D_Sav+D₊)` to `0.0`; and the
local margin `h_n − ⟨K₊⁺c,c⟩` is exhibited **negative** on a random (non-good) mode — `m_soft` is
**not** automatic. **[T]**

### 2.2 The inequality `m_soft ≥ δ` and its constants `[T]`

A positive margin `δ` needs **(★)**:
```
E[h_n]  −  E⟨K₊⁺c,c⟩  −  E|⟨φ₋|c⟩|²/(gB)  ≥  δ .
```
The three pieces and their constants:
- **bare stiffness** `E[h_n] = S_bare` — the coarse-mode diagonal Hessian (`~βQ_c + m²`), an **H1 input**;
- **positive-fibre Schur reduction** `E⟨K₊⁺c,c⟩ ≥ 0` — the *non-negative* cost of integrating out the
  `K₊⪰0` fibre; bounded by `E‖c‖²/γ₊` **only off the kernel** (`γ₊` = least positive eigenvalue of `K₊`);
- **rank-one defect** `D_Sav = E|⟨φ₋|c⟩|²/(gB) ≥ 0` — the negative-mode enhancement, `∝1/(gB)`.

**Sharp matrix-element form `[T]`.** By Cauchy–Schwarz in the fibre (`‖φ₋‖=1`) with the `+1`-helicity
projection,
```
D_Sav = E|⟨φ₋|P_{(+1)}H_cf v⟩|²/(gB) ≤ E‖P_{(+1)}H_cf v‖²/(gB) .
```
So **(★)** is implied by the single operator inequality on the coarse subspace
```
E‖P_{(+1)}H_cf v‖²  <  (gB)·E[h_n] − (gB)·δ .     ((★_op))
```
This is **the precise inequality**: the helicity-matched, one-octave coarse–fibre coupling energy must
be dominated by the negative-mode gap `gB` times the bare stiffness, with a `δ`-margin.

### 2.3 The refined rank-one bound and the `2/λ₄` correction `[T]/[V-num]`

Projecting onto the negative-mode amplitude (`V(φ)=−½(gB)φ²+¼λ₄φ⁴`) gives the toy reduction
`D_Sav ↔ g²Var(φ²)` (`SAVVIDY-DTERM` §2), so **(★_op)** is bounded by
```
D_Sav ≤ g² Var(φ²; gB, λ₄) ,   Var evaluated at the EXACT finite c=gB .
```
**Honest correction to the corpus.** `SAVVIDY-NEGMODE-EXACT` states `Var(φ²) ≤ 2/λ₄`. The exact
quadrature (witness (3)) shows `2/λ₄` is the **`c→∞` asymptote, not an upper bound** — there is a
**peak `2.21/λ₄` at `gB≈4`** (`+10%`):

| `gB` | 0.1 | 0.5 | 1 | 2 | 4 | 10 | 100 |
|---|---|---|---|---|---|---|---|
| `Var(φ²)` (`λ₄=1`) | 0.575 | 0.723 | 0.956 | 1.557 | **2.209** | 2.022 | 2.000 |
| `/(2/λ₄)` | 0.29 | 0.36 | 0.48 | 0.78 | **1.10** | 1.01 | 1.00 |

So the correct bookkeeping is `D_Sav ≤ g²Var(gB,λ₄)`, and at the **head** (`gB=O(1)`, small-`c`) the
exact `Var` is **`0.2–0.5×(2/λ₄)`**, i.e. **2–3× tighter** than the `2g²/λ₄` used previously. (The
`2/λ₄` "bound" is only valid deep in the condensate `gB≫1`; the head is precisely the small-`c`
regime where it is far too lossy.)

---

## 3 · The theorem attempt at `β=2.5`, its positivity region, and where it fails `[C]/[O]`

Take `g²(β)=4/β`, `β=2.5 ⇒ g²=1.6`, and the measured `S_bare = E[S] ≈ 1.40` (`[V-num]`, stable
`L4→L6`). Define the (toy-reduced) margin
```
m_toy(gB,λ₄) = S_bare − g² Var(gB,λ₄)  ≤  m_soft   (bound; lossy at the head).
```

### 3.1 Positivity region `[V-num]`

`m_toy>0` at `β=2.5` (witness (4)):

| `gB\λ₄` | 0.5 | 0.8 | 1.0 | 1.2 | **1.6** | 2.0 | 2.4 | 3.2 |
|---|---|---|---|---|---|---|---|---|
| 0.3 | m<0 | **+** | **+** | **+** | **+** | **+** | **+** | **+** |
| 0.5 | m<0 | m<0 | **+** | **+** | **+** | **+** | **+** | **+** |
| 1.0 | m<0 | m<0 | m<0 | **+** | **+** | **+** | **+** | **+** |
| 2.0 | m<0 | m<0 | m<0 | m<0 | **+** | **+** | **+** | **+** |
| 4.0 | m<0 | m<0 | m<0 | m<0 | m<0 | m<0 | **+** | **+** |

The **threshold** `λ₄,min(gB)` (thinnest quartic with `m_toy>0`): `0.76` (`gB=0.3`), `0.85` (`0.5`),
`1.07` (`1.0`), `1.51` (`2.0`), `1.90` (`3.0`), `2.22` (`4.0`). Natural `λ₄~g²=1.6` sits **above** the
threshold for `gB≲2` (margin **positive**, e.g. `gB=1,λ₄=1.6 ⇒ m_toy=0.55`) and on/below it for
`gB≳3`. So a **conditional positive margin exists** in the region `{λ₄ ≥ λ₄,min(gB)}`.

### 3.2 `β`-dependence `[C]/[O]`

`soft = g²Var(gB,λ₄)`; its `β→∞` limit is an **RG input** (witness (5)):
- **`λ₄` fixed O(1):** `soft = g²Var → 0`, margin `→ S_bare>0` — **safe** `[C]`.
- **natural `λ₄=g²`** (quartic runs with the coupling), `gB` fixed: the **asymptote** `soft→2g²/λ₄ = 2`
  **exceeds** `S=1.4` (the `2/λ₄` regime), **but** the **finite-`gB` exact** value at `β=2.5` is
  `0.85<1.4` (pass). So the pass **hinges on small `gB`**, and at large `β` under natural running the
  margin drifts toward the marginal asymptote (`m_toy(β=50)≈−0.79` in the `gB=1` toy).

### 3.3 Where it fails — the exact obstruction `[T]/[O]`

The conditional margin `m_toy>0` is a **bound** (`m_soft ≥ m_toy`), valid only if `g²Var` saturates
`D_Sav`. But **(★★)** is decisive and exact:

> **Ceiling theorem `[T]`.** `m_soft = G″ − |D₊| ≤ G″`. Since `D₊≤0` is **unbounded below** (BL gives
> only the *sign*), `m_soft ∈ (−∞, G″]`. Measured `G″ ≈ +0.049 = 3.5%·E[S]` ⇒ **`m_soft ≤ 0.035·S`**,
> possibly `≤ 0`. Witness (7b).

Consequences (the precise failure):
1. **No O(1) margin.** Any positive soft-mode margin is at most the measured `G″`, i.e. **thin
   `O(3%)·S`** — not the `m_toy=0.55` of the `gB=1` corner. The measured `D/S≈0.97` forces `gB` into
   the marginal regime `gB≳2.1` (natural `λ₄`), where `m_toy ≲ 0.05` (`≤0` for larger `gB`),
   consistent with `(★★)`.
2. **Wrong carrier.** The observed positivity `G″>0` is delivered by the **BL term** `|D₊|`
   (`G″=m_soft+|D₊|`), not by the rank-one soft mode. Proving `m_soft>0` is *sufficient* for `G″>0`
   but **strictly stronger** than needed, and **not supported by the data**.
3. **Lossy bound at the head.** At `k~√(gB)` the overlap is genuinely `O(1)` (§1.2), so the toy bound
   `D_Sav ≤ g²Var` (which assumes the isolated-cell amplitude) is **lossy by an O(1) factor**; the
   matrix-element bound `E‖P_{(+1)}H_cf v‖²/(gB)` is the honest one and is exactly `(★_op)`, unproved.
4. **The exact missing step.** `(★_op)` at the head, `E‖P_{(+1)}H_cf v‖² < gB·E[h_n] − gB·δ`, is (a) an
   **H1 input** through `E[h_n]=S_bare` and `gB` (magnetic scale / condensate), and (b) accumulated
   **IR-marginally** (`H2` §2.4, `SAVVIDY-MULTICELL-RG`: `Σ ~ C_RG Var(φ²) ln ln(k_UV/Λ)/ln(k_IR/Λ)`),
   decided at the **Landau pole `k_IR→Λ`** where one-loop tools die.

**Precise obstruction (the requested `[O]`).** Not a divergence, and not even the `O(1)` overlap per
se — it is the **magnitude** of the BL buy-in `|D₊|`: the ceiling `m_soft ≤ G″` shows the rank-one
soft-mode margin is pinned **below the (thin, `L`-shrinking) measured `G″`**, so a *positive
soft-mode* margin is at best `O(3%)·S` and the positivity is really BL-carried. Proving the *needed*
object needs a **lower bound on `|D₊|`** (the non-Gaussianity of the positive fibre), which BL does
not give. `[O]`.

### 3.4 The two exact statements pin `gB` and close the soft-mode route `[T]/[V]`

The bound `m_soft ≥ m_toy = S − g²Var` and the ceiling `m_soft ≤ G″` must coexist, forcing
`m_toy ≤ G″ ⇔ g²Var(gB,λ₄) ≥ S − G″ ≈ 0.965·S`. At natural `λ₄=1.6` (witness (7c)):
```
gB=1.0: g²Var=0.852 (<0.965 S)    gB=2.0: 1.293 (<)    gB=2.5: 1.542 (>=)
```
So consistency requires **`gB ≳ 2.1`** — precisely the regime where `m_toy ≲ G″ ≈ 0.05` (and `≤0` for
`gB` larger). The optimistic `gB=1` corner (`m_toy=0.55`) is **excluded by the data** (`D_Sav≈S`).
**The two exact statements together pin the head's negative-mode curvature to `gB≳2.1` and, *at that*
`gB`, cap the soft-mode margin at `O(3%)·S`.** This is the closed loop: the rank-one structure is
self-consistent only in the thin-margin regime, and there it cannot certify `δ` of order `S`.

---

## 4 · Ledger, acid test, reading

| item | statement | label |
|---|---|---|
| overlap lemma `⟨φ₋|c⟩ = ⟨H_cf φ₋|v⟩` | one coarse–fibre Hessian matrix element (err 0) | **[T]** |
| overlap `O(1)` at head `k~√(gB)`; suppressed off it | scale separation / one-octave / helicity | **[T]** (LLL model / structural) |
| deep-IR `\|overlap\|²∝k²→0`, UV `e^{−k²/gB}`, gauge sector `=0` | transversality / LLL / `+1`-helicity | **[T]** |
| margin identity `G″ = (S−D_Sav)+|D₊| = m_soft+|D₊|` | exact (witness (1)) | **[T]** |
| inequality `m_soft≥δ ⟺` (★); matrix form (★_op) | exact rewrite | **[T]** |
| refined bound `D_Sav ≤ g²Var(gB,λ₄)` | exact finite-`c` quadrature | **[T]/[V]** |
| `Var(φ²)≤2/λ₄` **false**; peak `2.21/λ₄` at `gB≈4` | exact quadrature; small-`c` head is `≤0.5·(2/λ₄)` | **[T]** (correction) |
| conditional margin `m_toy>0` region `{λ₄≥λ₄,min(gB)}` | `β=2.5`, `S_bare=1.40` | **[C]/[V]** |
| `β→∞` safety iff `λ₄/g²→∞`; natural running marginal | `soft=g²Var` asymptote `→2g²/λ₄` | **[C]/[O]** |
| **ceiling** `m_soft = G″−|D₊| ≤ G″`; measured `≤0.035·S` | the decisive obstruction | **[T]** |
| bound+ceiling pin `gB ≳ 2.1` at the head (natural `λ₄`) | `m_toy ≤ G″`; data `D_Sav≈S` | **[T]/[V]** |
| positivity is **BL-carried**, not rank-one-carried | `G″=m_soft+|D₊|`, data `D/S≈0.97` | **[T]/[V]** |
| uniform pointwise `m_soft>δ>0` | false: ceiling + measured thin/shrinking `G″` | **[NG]** |
| head = H1 carried by `|D₊|` magnitude / `(★_op)` | IR-marginal at Landau pole | **[O]** |
| **(★) at the head** | **not closed** | **[O]** |

**Acid test (free field), for every tool used.**
- **Overlap lemma / split:** `B=0` ⇒ no Savvidy background ⇒ no negative mode `φ₋` ⇒ `D_Sav=0`,
  `D₊=0` (Gaussian fibre) ⇒ `m_soft=S`, `G″=S`. At the IR mode `Q_c~4π²/L²→0`, `S→0`. The tool
  **manufactures no gap** on the free field. **Passes** (witness (6)).
- **Refined `Var` bound:** at `β→∞` with `λ₄` fixed, `g²Var→0`; the bound asserts no positivity on a
  free field (it bounds `D_Sav`, does not manufacture `G″>0`). **Passes.**
- **Ceiling `m_soft≤G″`:** on the free field `D₊=0` and `m_soft=G″=0`; equality, no false positivity.
  **Passes discriminatingly.**

**Reading.** The rank-one soft-mode structure is *fully* exploited: overlap exact (lemma) and O(1) at
the head; the margin identity exact; the finite-`c` bound refined (and the `2/λ₄` asymptote error
corrected). What it delivers is a **conditional** positive margin in a computable `(gB,λ₄)` region,
and — decisively — the exact ceiling `m_soft = G″ − |D₊| ≤ G″`. That ceiling is the honest verdict:
the rank-one soft mode **cannot** supply an O(1) positive margin, because the measured positivity is
carried by the **BL term `|D₊|`**, and `|D₊|`'s *magnitude* (not its sign) at the head is the exact
missing step `(★_op)` = H1. **Clay 0%. No QED.**

---

## 5 · Witness — `research/v3/w7/sm2_rank_one.py`

Deterministic, `py -3.12 -I`, numpy only. Quoted output (abridged):

```
(1) OVERLAP LEMMA:  <phi_-|c> = <H_cf phi_-|v>            err = 0.00e+00
    split <K^-1 c,c> = <K_+^+c,c> - <phi,c>^2/gB           err = 7.11e-15
    toy local margin h_n - <K_+^+c,c> = -7.79  (random mode: m_soft NOT automatic)
(2) OVERLAP vs k/sqrt(gB):  r=0.1 -> 0.0098 ; r=0.79 -> 0.2059 (PEAK, O(1)) ;
    r=2 -> 0.0147 ; r=3 -> 0.0001   (deep-IR -> 0, one-octave O(1), UV decay)
(3) EXACT Var(phi^2): c=0.5->0.723, c=1->0.956, c=4->2.209 (PEAK >2/lam), c=100->2.000
    => 2/lam4 is the c->inf ASYMPTOTE, not a bound (+10% peak); head (small c) 0.2-0.5 x it
(4) beta=2.5 g^2=1.6 S_bare~1.40:  m_toy = S - g^2 Var(gB,lam4)
    gB=1,lam4=1.6 -> 0.548 (YES) ; gB=2,lam4=1.6 -> 0.107 (YES) ; gB=2,lam4=0.8 -> -2.005 (NO)
    lam4,min(gB): 0.76(gB=.3) 1.07(1.0) 1.51(2.0) 2.22(4.0)
(5) beta-dep: lam4 fixed -> g^2Var->0 (SAFE) ; natural lam4=g^2 -> asymptote 2g^2/lam4=2 > S
(6) ACID TEST: B=0 -> D_Sav=0, m_soft=S=0, G''=0  (free massless; passes)
(7) MEASURED (beta=2.4,L=8,Ncfg=400):  D/S mean=0.9651 ; G'' mean=+0.0490 (3.5% S) ; 7/12 G''>0
(7b) CEILING: G''=m_soft+|D_+| => m_soft <= G'' = +0.0490  (positivity is BL-carried)
(7c) PIN: need g^2Var >= S-G'' ~ 0.965S: gB=1->0.85, gB=2->1.29, gB=2.5->1.54 => gB>~2.1
```

Readings: (1) the overlap lemma and split are exact; the local margin is sign-indefinite (not
automatic); (2) the overlap is O(1) at the head, suppressed off it by the three named mechanisms;
(3) the exact `Var` corrects the `2/λ₄` asymptote and is 2–3× tighter at the head; (4) a conditional
positive margin exists for natural-to-stiff `λ₄`, `gB≲2`; (5) `β→∞` safety needs `λ₄/g²→∞`; (6) acid
test clean; (7)/(7b) the measured positivity is thin and **BL-carried** — `m_soft ≤ G″ = 3.5%·S`.
**Not a proof of H1 — a witness for the exact overlap, the exact margin, the refined bound, and the
exact ceiling that decides the task.**

---

## 6 · References (internal, this session)

- `H2-DEFECT-BOUND.md` §0–§2 (rank-one split `D=D₊+D_Sav`, `D_Sav=E|⟨φ₋|c⟩|²/(gB)`, one-octave/deep-IR/
  `2g²/λ₄` routes, Fisher margin `I(ρ)`).
- `SAVVIDY-DTERM.md` §1–§3 (exact `D≤D_Sav`, toy `Γ″=κ−g²Var(φ²)`, deep-IR vs crossover).
- `SAVVIDY-NEGMODE-EXACT.md` §1–§4 (exact `Var(φ²)`, `±v` saddle, `c`-independence asymptote).
- `SAVVIDY-UNIFORM-SOFTENING.md` §2–§3 (one-octave `(L,a)`-uniform, marginal `O(1)` in `β`).
- `SAVVIDY-2LOOP.md` §2–§3 (sign protection, small 2-loop constant shift).
- `H1-RG-DESCENT.md` §4 (head = H1, soft-mode failure `η_j≥1`); `SAVVIDY-MULTICELL-RG.md` §2–4
  (IR-marginal log-log, Landau pole).
- Brascamp–Lieb (1976), Prékopa (1971); Landau level / Savvidy (standard).
- Witness `research/v3/w7/sm2_rank_one.py` (deterministic; output quoted §5).

*No third-party sources downloaded; all external inputs are standard and cited from the corpus.*
