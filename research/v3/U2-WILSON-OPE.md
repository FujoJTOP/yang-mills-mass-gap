# U2-WILSON-OPE — the Wilson-loop shape derivative and the two-`tr F²` OPE are **exact identities for the DIAGONAL worldsheet condensate / string tension**, and they give **no lower bound on the tube overlap `w₁`**; the atom is an *off-diagonal spectral-edge residue* the loop cannot see, and 2D YM (`σ>0`, `f₁²=0` exactly) plus the free field (`σ=0`, `f₁²=0`) are the decisive acid tests

> **2026-11-24 · Wave U2 — the Wilson-loop attack on the nonperturbative binding theorem.**
> The T wave left one object: a uniform lower bound on the atom fraction `w₁=f₁²/M₀` of the minimal four-wall flux tube (`T-WAVE §5`, `T2-TUBE-OVERLAP §8`). The tube *is* a flux tube; its overlap might be expressible through Wilson-loop geometry. This note writes the three candidate identities — (i) the **area/shape derivative** of `⟨W(C)⟩` as the insertion of `tr F²` on the worldsheet, (ii) the **short-distance OPE** of two `tr F²` insertions, (iii) the **trace-anomaly / conformal Ward** relation between `f₁`, the glueball mass `m_G`, and a decay constant — and decides, for each, whether it lower-bounds `w₁`, and in which direction positivity of the string tension acts.
> **Reads.** `T-WAVE-SYNTHESIS.md`, `T2-TUBE-OVERLAP.md`, `T5-CENTER-VORTEX.md`, `T6-BINDING-ADVERSARIAL.md`, `R5-EXACT-ANCHORS.md`; `unconditional/PROOF-P3-UNCONDITIONAL.md`; `docs/paper-v2/ch04-gap.tex`, `ch02-squeeze.tex` (`eq:squeeze`, `§schur`).
> **Witness.** `research/v3/w7/u2_wilson_ope.py` (`py -3.12 -I`, numpy only, deterministic, <1 s; output quoted §7). No third-party downloads.
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional/model · **[NG]** no-go for a named class · **[O]** exactly the missing step · **[V]** numerical witness here · **[std]** standard textbook tool cited from the corpus, not re-fetched. **No QED. Clay 0%.**

---

## 0 · Verdict (up front; no QED, Clay 0%)

1. **The identities are exact, and they are exact statements about the wrong object `[T]/[NG]`.** The area derivative of the Wilson loop is the insertion of `F_{μν}` on the contour; its second (coincident, contracted) derivative is the insertion of `tr F²` on the worldsheet (Mandelstam/Polyakov operator identity, modulo the contact/self-energy term that is exactly the UV renormalisation of `tr F²`). The OPE of two `tr F²` is the Wilson OPE whose coefficients are the **continuum** (unit operator `∼x⁻⁸`, `tr F²` operator, …). The trace anomaly gives an **equality** relating `f₁`, `m_G`, and a glueball decay constant `f_G`. All three are **exact** `[T]` (the 4D area law itself is `[C]/[T-standard]`). **None contains `w₁`.**

2. **The shape derivative and the OPE control the DIAGONAL (on-loop / bulk) object, not the OFF-DIAGONAL residue that `w₁` is.** The second area derivative gives `⟨tr F²(x) W(C)⟩_c` — the **local energy density of the string on its worldsheet**, i.e. (through the area coefficient) the **string tension `σ`**. The OPE gives `M₀ = ∫dν`, the **total spectral weight**, whose short-distance part is the perturbative continuum. But `w₁ = f₁²/M₀` with `f₁ = ⟨G|tr F²|Ω⟩` is the **vacuum-to-lightest-state matrix element** — an *off-diagonal* residue at the **bottom of the spectrum**. A diagonal condensate and a continuum are not the residue. **[T]**

3. **Direction: positivity of `σ` does NOT lower-bound `w₁`; the only bound-shaped content is an upper bound `[NG]`.** The OPE fixes the continuum and hence gives a *lower* bound on the denominator `M₀`; the trace anomaly is an equality, not an inequality. A lower bound on `M₀` plus no lower bound on `f₁²` can at most produce an **upper** bound on `w₁ = f₁²/M₀`. This is the exact reprise of `ch02 §schur` ("upper-bound machine, lower-bound target," family 14): every diagonal/moment/OPE functional is Schur-concave and bounds the atom from **above**, never below. **[NG]** for the class of shape-derivative/OPE-functionals.

4. **Decisive test — 2D pure YM is exactly solvable and shows `σ>0` with `f₁²=0` `[T]+[V]`.** In 1+1D the rep-`R` Wilson loop obeys an **exact area law** `⟨W_R⟩=e^{−σ_R A}` with `σ_R ∝ C₂(R)>0`; the area (shape) derivative is exactly `−∂_A log⟨W_R⟩=σ_R>0`. Yet `tr F²=tr E²` is the **Casimir Laplacian, diagonal in the character basis**, so its **on-loop** (diagonal) value is `C₂(R)/ℓ>0` (the string energy density, `∝σ_R`) while its **vacuum-to-atom** (off-diagonal) value is `⟨χ₀|tr F²|χ_R⟩=C₂(R)⟨χ₀|χ_R⟩=0` by character orthogonality. **The string tension is positive and the atom is exactly zero in the same solvable model.** This is the Wilson-loop coordinate form of T5's charged-vs-neutral blindness (`σ` in the charged sector, `tr F²` neutral). **[T]+[V]**

5. **Free field: `σ=0` and `f₁²=0` `[T]+[V]`.** No area law (perimeter/Coulomb), and `tr F²` of a free field is a two-particle cut (absolutely continuous, no delta) ⇒ `w₁=0`. The route passes the acid test — it manufactures no floor. **[T]**

6. **Whether the identity is exact or heuristic.** Split cleanly: the **operator identities** (area derivative ↔ `F`, second ↔ `tr F²`; OPE; trace anomaly) are **exact `[T]`**; the **4D area law** is `[C]/[T-standard]` (exact in 2D and at strong coupling, conjectural in the 4D continuum); the **step from these to `w₁`** is the **failing, heuristic step `[NG]`**. The identity "expresses `w₁` as a shape derivative" is *false*; it expresses the **diagonal condensate / string tension** as a shape derivative.

7. **The exact missing estimate is unchanged `[O]`:** a uniform `w₁ = f₁²/M₀ ≥ w₀>0` through `β=2.5`, equivalently `Σ_{n≥2}w_n ≤ 1−w₀` = isolation `I2` = `H*` (`T2 §8`). The Wilson-loop route does not supply it: to convert the **diagonal string data** into the **off-diagonal spectral-edge residue** one needs the lightest-glueball pole and its isolation — which *is* the gap. The route is therefore **circular, not a new species**; its only genuinely new object (the linearised Makeenko-Migdal loop equation = a loop-space spectral problem) has the glueball spectrum as its spectrum, so it renames `H*` rather than escaping it. **[O]/[C]**

**One line.** *The Wilson-loop shape derivative and the `tr F²` OPE are exact identities for the diagonal on-loop condensate and the string tension `σ`, but the tube overlap `w₁` is the off-diagonal vacuum-to-glueball residue at the bottom of the spectrum; positivity of `σ` does not lower-bound it — 2D YM has `σ>0` and `f₁²=0` exactly, the free field has `σ=0` and `f₁²=0` — and the only bound-shaped output of the route is an upper bound on `w₁` through a lower bound on `M₀`; the missing estimate is still `w₁≥w₀` = isolation = `H*`. No QED. Clay 0%.*

---

## 1 · The shape-derivative identity: `tr F²` as an area derivative of `W(C)`

### 1.1 The area derivative is the `F`-insertion `[T]`

Let `W(C)=N⁻¹tr P exp(∮_C A)` in representation `R`, and let `σ_{μν}(x)` be the worldsheet
area element (an infinitesimal deformation of the contour at `x` in the `μν` plane). The
**Mandelstam–Polyakov area derivative** is the exact operator identity
```
   δ/δσ_{μν}(x) W(C)  =  i g N⁻¹ tr[ F_{μν}(x) P exp_{C}(A) ] .                    (1.1)  [T]
```
The right side is the field strength at the point `x` of the worldsheet, path-ordered into the
loop. Equation (1.1) is the statement that *deforming the loop in the `μν` plane inserts `F_{μν}`
there*.

### 1.2 The second area derivative is the `tr F²` insertion `[T]` (mod contact)

Applying (1.1) twice and taking the coincidence limit `x→y` with the two area elements
contracted (`μν` with `μν`),
```
   (δ²/δσ²)_{μν,μν}(x) ⟨W(C)⟩  =  −g² ⟨ tr F²(x) W(C) ⟩  +  (contact) ,            (1.2)  [T] (mod contact)
```
where the **contact** term is the short-distance self-energy of the insertion on the loop — the
delta supported on `C` coming from `∂_μ F_{μν}`. This contact is *exactly* the UV divergence of
`⟨tr F²⟩` identified in T4 (`χ∼a⁻⁴`, needs renormalisation at a fixed physical scale). So (1.2)
is exact as an **operator** identity, with the **finite renormalised** statement `[C]`.

> **What (1.2) is, and what it is not.** The object it produces is
> `⟨tr F²(x) W(C)⟩_c / ⟨W(C)⟩` — the *conditional* expectation of `tr F²` **on a fixed loop**, i.e.
> the **local energy density of the string on its worldsheet**. In the string picture this is
> `≈σ` pointwise, and its integral over the cross-section is the tension. It is a **diagonal**
> (on-loop, normalised by the same loop) object. The atom fraction `w₁=f₁²/M₀` is a
> **vacuum-to-state matrix element** `⟨G|tr F²|Ω⟩` divided by a **vacuum two-point function** — an
> **off-diagonal, unnormalised** object. The normalisation by `⟨W(C)⟩` in (1.2) is precisely what
> makes the on-loop insertion nonzero while the vacuum-to-atom matrix element can be (and in 2D *is*)
> zero. **The shape derivative cannot produce the off-diagonal residue.**

### 1.3 The Makeenko–Migdal loop equation — a closed equation for `⟨W⟩`, not for `w₁` `[C]`

The loop equation is the exact Schwinger–Dyson identity obtained by letting `∂_μ` in (1.1) hit the
field strength, `D_μF_{μν}=J_ν` (the loop current),
```
   ∂_μ^x (δ/δσ_{μν}(x)) ⟨W(C)⟩  =  g² ∮_C dy_ν δ⁴(x−y) ⟨W(C_{xy})⟩⟨W(C_{yx})⟩ .          (1.3)  [std]
```
It is a **closed non-linear equation for `⟨W⟩`** (Makeenko–Migdal, in the large-`N` planar limit). Its
content for `tr F²` is again the *contracted coincidence* of §1.2 — the **diagonal** energy
density — and the loop-splitting RHS. It is *not* an equation for the spectral residue `f₁²`.
In 4D the equation is **not rigorously solved** (it is exactly as hard as the master field / the
gap); so even as a route to the spectrum it is `[C]` at best. See §5 for the one structurally new
object it offers (the loop-space linearisation) and why it does not help.

---

## 2 · The short-distance OPE of two `tr F²` insertions

### 2.1 The Wilson OPE and its coefficients `[T]`

Two `tr F²` (dimension-4) operators have a short-distance Wilson OPE
```
   tr F²(x) tr F²(0)  =  C_1(x) 1  +  C_{F²}(x) tr F²(0)  +  …  ,                        (2.1)
```
- `C_1(x) ∼ c₀ / x⁸` — the **unit-operator** coefficient; positive, fixed, this is the perturbative
  **two-particle continuum** (the free-field `tr F²` two-point function).
- `C_{F²}(x) ∼ (anomalous dimension) · log(x²μ²) / x⁴` — the coefficient of `tr F²` itself,
  controlled by the trace-anomaly anomalous dimension.
The vacuum matrix element of (2.1) is the two-point function `⟨tr F²(x) tr F²(0)⟩`, whose spectral
representation is `∫dν(s)` with `M₀=∫dν`, `M₁=∫s dν`, etc.

### 2.2 The OPE fixes the CONTINUUM, not the atom `[T]`

The coefficients `C_1, C_{F²}, …` are **short-distance** data: they fix the **UV / continuum** part of
`ν`. The atom `f₁²=ν({m₁²})` is a **long-distance spectral-edge residue** — the amplitude for the
lightest *state* to be created by the local operator. It is **not one of the OPE coefficients**
(a state is not a local operator in the OPE list, except through its interpolating field, whose
coefficient is not `f₁²`). Concretely,
```
   M₀ = ∫dν = f₁² + ∫_{cont} dν  =  f₁² + M₀^{OPE-continuum} + (higher-operator corrections) ,
```
so the OPE fixes `M₀^{continuum}` and hence gives `M₀ ≥ M₀^{cont} > 0` (a **lower** bound on the
denominator), while the atom `f₁²` sits *above* the OPE's reach. This is exactly the
`ch02 §schur` statement that finite moment/OPE data is blind to the spectral edge.

---

## 3 · The string tension: what it bounds, and in which direction

### 3.1 The area law and the diagonal condensate `[C]/[T-standard]`

The area law `⟨W_R(C)⟩ = exp(−σ_R·Area − μ·Perimeter + …)` defines the string tension `σ_R`. Its
area (shape) derivative is (1.2) at long distance,
```
   −∂_{Area} log⟨W_R(C)⟩ = σ_R + (curvature/perimeter corrections → 0 for large loops) .       (3.1)
```
So `σ_R` **is** the diagonal `⟨tr F²⟩` per unit area of the worldsheet, renormalised. The trace
anomaly ties the *trace* of the stress tensor to `tr F²`,
```
   T^μ_μ = (β(g)/2g) tr F² ,   ⟨0|T^μ_μ|G(p)⟩ = m_G² ⟨0|φ_G(p)⟩  ⟹  f₁ = (2g/β) m_G² f_G ,   (3.2)  [T]
```
where `f_G=⟨0|φ_G|G⟩` is the glueball decay constant of the interpolating field `φ_G`. Equation
(3.2) is an **equality relating three numbers**, not an inequality.

### 3.2 The direction theorem: no lower bound on `w₁` `[NG]`

Collecting §2.2 and §3.1:
- **OPE** gives `M₀ ≥ M₀^{cont}>0` — a **lower bound on the denominator** of `w₁=f₁²/M₀`.
- **Trace anomaly** (3.2) is an **equality**: `f₁ ∝ m_G² f_G`, no sign forced.
- **Area law** gives `σ_R>0` — a **diagonal** energy density, not the residue.

> **Direction theorem (the Wilson-loop route is upper-bound-shaped). `[NG]`**
> *The shape derivative and the OPE bound `σ`, the diagonal on-loop `⟨tr F²⟩`, and the continuum
> `M₀^{cont}` from below, and relate `f₁` to `m_G, f_G` by an equality. They impose **no lower bound
> on `f₁²>0`**. Since `M₀` is bounded from below and `f₁²` is not, the route can at most yield an
> **upper** bound `w₁ ≤ f₁²/M₀^{cont}` (and only if `f₁²` is independently bounded above). A
> lower-bound target cannot be reached by an upper-bound machine (`ch02 §schur`).*

**This is the answer to the task's direction question:** positivity of the string tension gives
**neither a clean lower nor a clean upper bound** on `w₁`; what it genuinely controls is the
**diagonal** condensate/`σ` and the **denominator** `M₀`. In bound-shaped terms it is an
**upper**-bound ingredient, never a lower one.

---

## 4 · Tests: 2D YM (solvable) and the free field

### 4.1 2D pure YM — exact area law with `σ>0`, atom exactly zero `[T]+[V]`

In 1+1D, `tr F²=tr E²`; the rep-`R` Wilson loop couples only to the character state `χ_R`
(`⟨χ_{R'}|W_R|χ₀⟩=δ_{RR'}`), so the area law is **exact** (Migdal/Kazakov):
```
   ⟨W_R(C)⟩ = e^{−σ_R A} ,    σ_R ∝ C₂(R) > 0   (R≠0) ,    −∂_A log⟨W_R⟩ = σ_R .           (4.1)  [T]
```
And `tr F² ∝ (1/ℓ)(−Δ_G)` is the **Casimir Laplacian, diagonal in the character basis** (the
precise positive prefactor is a convention; only `∝C₂(R)` and its positivity are used):
```
   on-loop (diagonal):   ⟨tr F²⟩_{W_R}/⟨W_R⟩ = c·C₂(R)/ℓ > 0      ( energy density ∝ σ_R )     (4.2a)
   vacuum-to-atom:       ⟨χ₀|tr F²|χ_R⟩ = c·(C₂(R)/ℓ)·⟨χ₀|χ_R⟩ = 0  (R≠0)                     (4.2b)  [T]
```
some `c>0`, because `⟨χ₀|χ_R⟩=δ_{0R}` (Haar orthogonality; witness (C)). Therefore, **exactly**,
```
   σ_R > 0    AND    f₁² = 0    ( w₁ = 0 )       in 2D pure YM .                               [T]+[V]
```
The area (shape) derivative is nonzero and positive; the atom is zero. **The 2D solvable model is the
decisive counterexample to "`σ>0` ⇒ `w₁>0`"** — and it is the Wilson-loop form of T5 (`σ>0` in the
*charged* sector, `tr F²` *neutral* and blind). Note also the **diagonal vacuum condensate
`⟨tr F²⟩=C₂(0)=0`** here: in 2D, `σ` and `tr F²` are almost decoupled — the string tension lives in
the charged sector and `tr F²` sees zero in the neutral vacuum. (This is why `R5 §1.1` found the 2D
`tr F²` atom zero while the 2D `0⁺⁺` "glueball" is a *holonomy* state.)

### 4.2 Free field — `σ=0`, atom zero (acid test) `[T]+[V]`

A free gauge field has no area law (`σ=0`, perimeter/Coulomb), so (3.1) returns zero; and `tr F²` of a
free field is a **two-particle bilinear**, whose KL measure is the continuum `[4m²,∞)` (or `[0,∞)`)
with **no delta**, so `f₁²=ν({m²})=0` (`R5 §1.6`). The route returns `w₁=0`: **it manufactures no
floor**, which is the acid test.

### 4.3 3D / 4D — no new mechanism beyond T2 `[T-standard]/[C]`

The area law is `[T-standard]` in 3D and lattice 4D, `[C]` in the continuum. Nothing in the
shape-derivative/OPE structure changes with dimension: it still produces the diagonal `σ` and the
continuum `M₀`, still not the residue. So 3D/4D inherit the 2D no-go: the area law (`σ>0`) is
compatible with `f₁²` being anything from zero up.

---

## 5 · Exact or heuristic? And is there a new species? (honest)

**Exactness, item by item.**

| identity | status |
|---|---|
| area derivative `=` `F`-insertion (1.1) | **[T]** (Mandelstam/Polyakov operator identity) |
| second area derivative `=` `tr F²`-insertion (1.2) | **[T]** mod contact; finite renormalised form **[C]** (T4 UV) |
| Makeenko–Migdal loop equation (1.3) | **[std]** as an identity; **unsolved** in 4D `[C]` |
| Wilson OPE of two `tr F²` (2.1) | **[T]** (standard OPE); the atom is *not* a coefficient |
| trace anomaly `f₁=(2g/β)m_G²f_G` (3.2) | **[T]** as an equality (conformal Ward) |
| 4D area law `⟨W⟩=e^{−σA}` | **[C]/[T-standard]** (exact 2D/strong coupling; conjectural continuum) |
| **step to `w₁`** | **fails `[NG]`** — no lower bound |

So: **the operator identities are exact; the bridge to the atom is the heuristic/failing step.**

**Is there a new species? `[O]/[C]`** The charter (§1) asks for a genuinely different object, not an
old inequality re-dressed. Here the honest answer is **no new species**, and one partial exception:
- The shape-derivative + OPE route is the **Wilson-loop coordinate form** of two already-known
  obstructions: T5's charged-vs-neutral blindness (`σ` charged, `tr F²` neutral; 2D witness) and
  `ch02 §schur`'s diagonal-vs-edge (Schur-concave / upper-bound) split. It does not escape either.
- The **one** structurally new object is the **linearisation of the Makeenko–Migdal loop equation
  about the confining solution** — a **loop-space spectral problem** whose spectrum is the
  glueball spectrum and whose propagator's residue is a loop-space `f₁²`. But (a) the confining
  solution of (1.3) is unproven in 4D, and (b) its spectrum *is* the gap, so this is a **renaming of
  `H*`, not a weakening of it**. `[C]` on a hypothetical solution, `[O]` otherwise.

---

## 6 · The exact missing estimate (unchanged), and why this route is circular

**Missing estimate (stated once, exactly).**
> There exist `w₀>0` and `δ>0` such that, **uniformly in `(L,a)` and for all `β` through the physical
> `β=2.5`**, the normalised flux-tube overlap obeys
> ```
>   w₁(β;L,a) = f₁²/M₀ ≥ w₀ > 0 ,     equivalently  Σ_{n≥2} w_n ≤ 1 − w₀ ,
> ```
> with the **isolation margin** `Λ−m₁²≥δ` (`I2`) that makes the second level usable in the
> residue–level-gap lemma. By `T2 §8` this is `I2`, which is `H*`, which is the gap.

**Why the Wilson-loop route is circular (does not supply the estimate).** The chain is
```
   σ>0  (diagonal, on-loop)   ⇏   f₁²>0  (off-diagonal, edge)                            [NG]
   f₁²>0  =  ⟨G|tr F²|Ω⟩  needs  the lightest glueball pole  +  its isolation I2         [O]
   isolation I2  =  the gap  =  H*                                                       [O]
```
The conversion of the *diagonal string data* (`σ`, the on-loop `⟨tr F²⟩`, the OPE continuum) into the
*off-diagonal spectral-edge residue* `f₁²` **requires knowing the lightest-glueball pole and that it
is isolated** — which is precisely the gap one rung up. The trace anomaly (3.2) makes this explicit:
it *relates* `f₁` to `m_G` and `f_G`, but is an equality that `f_G=0` satisfies for any `m_G>0`
(the 2D witness: `f_G=0` while `m_G, σ_R>0`). So the Wilson-loop route reproduces the circularity
`R5 §3.4` / `E2-REGRESS` rather than breaking it.

**What *would* be needed (the precise gap).** A proof that the diagonal string data **determines the
off-diagonal edge with a definite sign and floor** — equivalently, that the lightest `0⁺⁺` pole of
`tr F²` is isolated and has residue `≥w₀`. No shape derivative, OPE coefficient, area law, or
trace-anomaly identity does this: they are all **diagonal / equality / continuum** statements, and
the target is an **off-diagonal / inequality / edge** statement.

---

## 7 · Witness — `research/v3/w7/u2_wilson_ope.py`

`py -3.12 -I u2_wilson_ope.py` (exit 0; numpy 2.4.6; deterministic, <1 s). Quoted output (abridged):

```
(C) SU(2) Haar character overlaps  <chi_0|chi_R> = delta_{0R}:
      <chi_0|chi_0> = +1.000e+00 ;  <chi_0|chi_1..3> = -9.4e-16, -1.5e-15, -2.1e-15   (orthogonality)
(A) 2D YM exact area law <W_R> = exp(-sigma_R A), sigma_R=(g^2/2)C2(R):
      fund C2=0.75 -> sigma=0.375>0 ; adjoint C2=2 -> sigma=1.0>0 ;
      -d/dA log<W_fund> = +0.375000 (constant) == sigma_fund : True
(B) 2D tr F^2 diagonal in characters -> on-loop POSITIVE, atom ZERO:
      fund: on-loop <trF^2>=C2/l=0.7500  ;  atom <0|trF^2|R> = 0.0000
      adj : on-loop = 2.0000              ;  atom = 0.0000
      ==> sigma_R>0 AND f1^2=0 coexist EXACTLY  (no-go; <chi_0|chi_R>=delta_0R is the mechanism)
(D) OPE fixes continuum: M_cont=1.8856, M0=f1^2+M_cont>=M_cont ; w1=f1^2/M0:
      f1^2=1.0 -> w1=0.3465 ; 1e-2 -> 0.0053 ; 1e-4 -> 5.3e-5 ; 0.0 -> 0.0000
      ==> w1 -> 0 at fixed continuum; OPE bounds the DENOMINATOR, cannot force a bottom atom
(E) free field: sigma=0 ; tr F^2 two-particle cut, no delta => f1^2=0  (ACID TEST PASSES)
(F) trace anomaly f1=(2g/beta) m_G^2 f_G is a RELATION; f_G=0 => f1=0 (2D: f_G=0 while m_G,sigma>0)
```

**Readings.** (C) the orthogonality `δ_{0R}` is the mechanism that kills the atom; (A) the shape
derivative is exactly `σ_R>0`; (B) the same solvable model has **positive string tension and zero
atom** — the decisive no-go; (D) the OPE's knowledge of the continuum cannot prevent `w1→0`; (E) the
free field returns zero (acid test); (F) the trace anomaly relates but does not bound. **Witness,
not proof; the theorems are §1–§6, the numerics illustrate the diagonal-vs-edge split and the exact
2D values.** *(The 2D `(g²/2)` prefactor is a convention; only `σ_R∝C₂(R)>0` is used.)*

---

## 8 · Ledger

| item | statement | label |
|---|---|---|
| area derivative = `F`-insertion | `δ/δσ_{μν}(x)W = ig tr[F_{μν}(x)P]` | **[T]** (Mandelstam/Polyakov) |
| second area derivative = `tr F²`-insertion | `−g²⟨tr F² W⟩_c` + contact; finite form = diagonal on-loop `⟨tr F²⟩` | **[T]** mod contact; **[C]** renormalised |
| Makeenko–Migdal loop equation | closed eqn for `⟨W⟩`; unsolved in 4D | **[std]/[C]** |
| Wilson OPE of `tr F²` | `C_1∼x⁻⁸` (continuum), `C_{F²}` (anomalous dim); atom is **not** a coefficient | **[T]** |
| OPE output | `M₀ ≥ M₀^{cont}>0` — lower bound on the **denominator** | **[T]** |
| trace anomaly | `f₁=(2g/β)m_G²f_G` — **equality**, `f_G=0` allowed for any `m_G>0` | **[T]** |
| string tension / area law | `σ_R>0` = **diagonal** on-loop energy density; 4D area law conjectural | **[C]/[T-standard]** |
| direction | `σ>0` controls `σ` & `M₀`; **no lower bound on `f₁²`**; at most an upper bound on `w₁` | **[NG]** |
| 2D YM exact | `⟨W_R⟩=e^{−σ_RA}`, `σ_R∝C₂(R)>0`; `tr F²` diagonal ⇒ on-loop `C₂/ℓ>0`, **atom 0** | **[T]+[V]** |
| free field | `σ=0`; `tr F²` two-particle cut ⇒ `f₁²=0` (acid test passes) | **[T]+[V]** |
| exact vs heuristic | operator identities exact; area law `[C]`; **bridge to `w₁` fails `[NG]`** | — |
| new species? | no — Wilson-loop form of T5 charged/neutral + `ch02 §schur` diagonal/edge; loop-space spectral problem renames `H*` | **[O]/[C]** |
| **uniform overlap floor** | `w₁≥w₀>0`, `(L,a)`-uniform, `β` through `2.5` — **not supplied by this route** | **[O]** |
| **missing estimate =?** | `w₁≥w₀` ⟺ `Σ_{n≥2}w_n≤1−w₀` ⟺ isolation `I2` ⟺ **`H*` = the gap** | **[O]** |
| **Clay** | unbuilt | **[O] · 0%** |

---

## 9 · Sharpest line

**The Wilson-loop OPE and the area/shape-derivative identity are exact, and they are exact
statements about the wrong object. The area derivative inserts `F_{μν}` on the contour; its
contracted coincidence is the `tr F²` insertion on the worldsheet — but normalised by the loop, so
it yields the *diagonal* on-loop energy density of the string, i.e. the string tension `σ`, and the
two-`tr F²` OPE yields the *continuum* `M₀`, not the atom. The tube overlap `w₁=f₁²/M₀` is the
*off-diagonal vacuum-to-lightest-glueball residue at the bottom of the spectrum*, which neither the
diagonal shape derivative nor the UV OPE can see; the trace anomaly relates `f₁` to `m_G` and a decay
constant `f_G` by an *equality* that `f_G=0` satisfies for any `m_G>0`. So positivity of the string
tension does **not** lower-bound `w₁`: it controls `σ` and the denominator `M₀` (a lower bound), which
can only produce an **upper** bound on `w₁` — the `ch02 §schur` "upper-bound machine, lower-bound
target." The solvable 2D model settles it exactly: the area law gives `σ_R>0` while `tr F²`, being
the diagonal Casimir Laplacian, has on-loop value `C₂(R)/ℓ>0` and *vacuum-to-atom* value
`C₂(R)⟨χ₀|χ_R⟩=0` — positive string tension and zero atom in the same theory (the Wilson-loop form of
T5). The free field gives `σ=0` with `f₁²=0`, passing the acid test. The identity is exact; the
bridge from `σ` to `w₁` is the failing step. The exact missing estimate is unchanged —
`w₁≥w₀` uniformly through `β=2.5`, equivalently `Σ_{n≥2}w_n≤1−w₀` = isolation `I2` = `H*` = the gap —
and the Wilson-loop route reproduces the circularity rather than breaking it. **No QED. Clay 0%.**

---

*Refs.* `T-WAVE-SYNTHESIS.md` (§5 next object), `T2-TUBE-OVERLAP.md` (§8 missing estimate), `T5-CENTER-VORTEX.md` (charged-vs-neutral, 2D `σ>0`, `f₁²=0`), `T6-BINDING-ADVERSARIAL.md`, `R5-EXACT-ANCHORS.md` (§1.1 2D atom zero, §1.6 free field, §3.4 residue–level-gap circularity); `unconditional/PROOF-P3-UNCONDITIONAL.md`; `docs/paper-v2/ch02-squeeze.tex` (`eq:squeeze`, `§schur` diagonal-vs-edge), `ch04-gap.tex` (`sec:savvidy`, gap anatomy). **Standard (textbook, cited from corpus, not re-fetched):** Mandelstam, Phys. Rev. **175** (1968) 1580 / Polyakov, Nucl. Phys. B **164** (1980) 171 (area derivative); Makeenko–Migdal, Phys. Lett. B **88** (1979) 135 / Nucl. Phys. B **188** (1981) 269 (loop equation); Wilson, Phys. Rev. **179** (1969) 1499 (OPE); Adler–Bardeen (trace anomaly) / Collins–Duncan–Joglekar, Phys. Rev. D **16** (1977) 438; Migdal, Sov. Phys. JETP **42** (1976) 413 / Kazakov–Kostov, Nucl. Phys. B **176** (1980) 199 (2D YM exact Wilson loop); Källén (1952)/Lehmann (1954). **Witness:** `research/v3/w7/u2_wilson_ope.py` (deterministic, numpy only, exit 0). **No third-party downloads; no `[C]` upgraded to `[T]`; no QED; Clay 0%.**
