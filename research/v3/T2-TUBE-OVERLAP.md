# T2-TUBE-OVERLAP — the minimal connected `¼tr F²` polymer is a four-walled flux tube; its overlap with the lightest `0⁺⁺` is a **normalised projection fraction** `w₁=f₁²/M₀ ∈[0,1]` that **does not collapse with tube length** but whose **certification does**; RP and PF cannot lower-bound it; the only proven floor is strong coupling, and at `β=2.5` the cluster expansion is `uK≫1` outside its radius; the exact missing estimate is a uniform floor `w₁≥w₀`, which by L2 **is** isolation/`H*`

> **2026-11-24 · Wave T2 — the binding criterion (the next species named by the R wave, `R-WAVE §7`).**
> R5 left one object: a **rigorous flux-tube/polymer lower bound** on the minimal connected `tr F²`
> overlap at `β=2.5`. This note writes that polymer and its overlap **explicitly**, asks what
> reflection positivity / Perron–Frobenius / the transfer matrix can and cannot give, and decides
> whether a positive bound is **uniform in the volume** and **survives at `β=2.5`**, or whether it
> **collapses with the tube length**.
> **Reads.** `R-WAVE-SYNTHESIS.md`, `R5-EXACT-ANCHORS.md`, `R6-ATOM-ADVERSARIAL.md`;
> `unconditional/PROOF-P3-UNCONDITIONAL.md`, `P3-FOURWALL-DERIVATION.md`, `P3-DRESSED-MASS.md`,
> `PROOF-GAUSS-INVARIANT-PF.md`, `G2D-ISOLATION-I2.md`, `GAP-G2-edge.md`, `BS-BLIND.md`,
> `STATE-DENSITY-I2.md`, `E2-REGRESS-CIRCULARITY.md`; `v3/SM2-RANK-ONE.md`, `H2-DEFECT-BOUND.md`.
> **Witness.** `research/v3/w7/t2_tube_overlap.py` (`py -3.12 -I`, numpy only, deterministic, ~0.2 s;
> output quoted §6). No third-party downloads.
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional/model · **[NG]** no-go for a
> named class · **[O]** exact missing step · **[V]** numerical witness here · **[std]** standard tool
> used as textbook. **No QED. Clay 0%.**

---

## 0 · Verdict (up front; no QED, Clay 0%)

1. **The minimal polymer and its overlap are explicit `[T]`.** The minimal connected, gauge-invariant,
   nonvanishing `¼tr F²` polymer linking the insertions at times `0,t` is the **straight four-walled
   tube** `B_T`: `4T` fundamental temporal plaquettes (`T=t/a`), cross-section the 4-link boundary of
   the source spatial plaquette, weight `w(B_T)=a₀u^{4T}`, `u=c_f/N`, `a₀=(2N)^{-p}s>0` (`P3-FOURWALL
   §1–§2`). Its **overlap** with the lightest `0⁺⁺` is `f₁=⟨G|O|Ω⟩=⟨G|φ⟩`,
   `f₁²=a₀(1−O(u))` **at strong coupling** — a **modulus-square of a positive `f⊗f̄→1` recoupling.**
   The honest, scheme-free form of "the overlap" is the **normalised atom fraction**
   ```
     w₁ := f₁²/M₀ ,   M₀ := Var(O) = ⟨O²⟩_c = ‖φ‖² ,   φ := Q\,O|Ω⟩ ,
   ```
   a **dimensionless number in `[0,1]`** (the fraction of the source's spectral weight in the isolated
   tube). It is **invariant under the multiplicative transport** `O→Z_O O` (numerator and `M₀` both
   scale by `Z_O²`) — this **sharpens R6-c**: the scheme-free floor is `w₁≥w₀>0`, not the
   `Z_O`-covariant `f₁²≥r₀`. **[T]**

2. **RP, PF and the transfer matrix do NOT give the lower bound `[NG]`.** Reflection positivity gives
   `dν≥0` and log-convexity of `C(t)` (⇒ `m_eff(t)↓m₁` *from above*), and Perron–Frobenius gives a
   strictly-positive vacuum and a Krein–Rutman *gap* — but **neither touches `w₁`**: log-convex `C(t)`
   is fully compatible with `w₁=0` (a pure cut), and PF is **vacuum-anchored** while the connected
   source `φ` is sign-indefinite, so the edge `|G⟩` lives in `Ω^⊥` where the one RP cone degenerates
   (`PROOF-GAUSS-INVARIANT-PF`, `BS-BLIND`, `G2D-ISOLATION-I2`). Both are **acid-test-blind**: the free
   massless field satisfies every consequence and has `w₁=0`. **No positivity/transfer-matrix functional
   lower-bounds the overlap — at any `β`, including `2.5`.** **[NG]** for the named class.

3. **"Collapse with tube length" — the object does not collapse; the *certification* does `[T]`.** The
   residue `w₁` is a property of the `t→∞` spectrum, hence **`T`-independent** on the cylinder; the
   strong-coupling prefactor `a₀` is `T`-independent and the dressed defect `D_fin(T)→1`
   (`P3-DRESSED-MASS`). At most roughening (transverse wandering) softens it **logarithmically** in 4D
   and is resummed into the **mass**, not the residue `[C]`. **But** the finite-`T` estimator
   ```
     ŵ₁(T) := e^{E₁ t}\,\bar C(t) = w₁ + Σ_{n≥2} w_n\,e^{−(E_n−E₁)t}   is DECREASING in t and ≥ w₁.
   ```
   So **every finite tube length OVER-estimates the overlap**; the true value is the limit, never
   certified below from finite `T` (and, by `R3` Stone–Weierstrass, never from finite `L`). The
   collapse is in **certifiability**, not in the object. **[T]** (witness (B)).

4. **At `β=2.5` there is NO proven positive uniform bound — and the only rigorous method is
   inapplicable there `[T]/[O]`.** The single unconditional floor is P3's strong-coupling residue
   (`f₁²≥r₀(β)>0`, `L`-uniform, `β<1/48`; `P3-DRESSED-MASS` §4). Its method is the Kotecký–Preiss
   cluster expansion, whose convergence needs `u·K<1` with `K≤CN`; at `β=2.5`, `u=I₂/I₁=0.719` and
   `uK∼10²–10³≫1`, **far outside the radius** (witness (F)). So the floor is not merely *untransported*
   — it is **not even defined** at `β=2.5`. The **object** (existence of a positive atom) is protected
   by confinement + the trace anomaly (`R6`), and the physical scaling `f₁²∼Λ⁶` plus the `β=2.2` data
   (`f₁²_bare=7.14`, `R5`) indicate it does **not** physically collapse `[V]/[C]` — but **no uniform
   positive lower bound is proven.** **[O]**

5. **The exact missing estimate + falsification test.** A **uniform (through `β=2.5`, in `(L,a)`)
   lower bound `w₁=f₁²/M₀ ≥ w₀>0`** — equivalently an **upper bound on the continuum leakage**
   `Σ_{n≥2}w_n ≤ 1−w₀`, i.e. *the source is mostly one tube.* By the residue–level-gap lemma L2 this is
   **identical** to a lower bound on the second level `Λ=E₂` (`w₁≥(Λ−Ē)/(Λ−m₁²)`, `Ē=M₁/M₀`), i.e. the
   **isolation `I₂`** — the `O(u)`-decoration sign/amplitude continued past the KP radius — i.e. **`H*`**
   = the gap. **Falsification test:** measure `w₁(L,T)` by GEVP + two-state fit of the zero-momentum
   `0⁺⁺` `tr F²` correlator; the floor holds iff `w₁` reaches a **plateau** that is `T`-independent and
   `L`-independent with value `>0`; it is **falsified** by a `T`-drift (`∂_T w₁<0` persisting) or an
   `L`-drift (`∂_L w₁<0`), i.e. leakage into the two-tube cut. **(Predicts: plateau, `w₁=O(0.2)–O(0.5)`,
   from the `β=2.2` pilot.)** **[O]/[V-plan]**

**One line.** *The flux-tube overlap is a `[0,1]` projection fraction `w₁=f₁²/M₀`; it does not collapse
with tube length (the object is `T`-independent), but no finite `(T,L)` certifies it and RP/PF cannot
lower-bound it; the only proven floor lives at `β<1/48` inside the KP radius, and `β=2.5` is `uK≫1`
outside it; the exact missing estimate is a uniform `w₁≥w₀`, which is `I₂`, which is `H*`, which is the
gap. No QED. Clay 0%.*

---

## 1 · The minimal connected `tr F²` polymer and its overlap, written out

### 1.1 The polymer `[T]`

Lattice: 4D pure SU(N) Wilson YM on the cylinder `T³_L×_aℤ`, `S=βΣ_p(1−N⁻¹Re tr U_p)`, character
expansion `e^{(β/N)Re tr U_p}=Σ_r c_r(β)χ_r(U_p)`, `c_f(β)=u+O(β³)`, `u:=c_f/N` (`=β/2N²` generic,
`β/N²` SU(2)-special; `P3-FOURWALL §4`). Carrier `O` = A1g zero-momentum single-**spatial**-plaquette
`¼tr F²`; `C(t):=⟨O(0)O(t)⟩_c`, `T=t/a`.

**Minimal polymer (proved, not asserted; `P3-FOURWALL §1`).** A connected gauge-invariant
nonvanishing polymer linking the two insertions has `|P|≥4T`:
- the source plaquette's **four** boundary links each carry an open fundamental pair, and `∫dU\,U=0`
  if uncovered, so `≥4` partner plaquettes (`≥4` flux strands);
- only **temporal** plaquettes advance time, and each of the `T` time-boundaries needs one per strand;
- merging strands into higher reps is neutral-to-worse in the cost gradng.

The **unique minimisers** `|P|=4T` are the **straight four-walled tubes** `B_T`: four side strips of
height `T`, one fundamental temporal plaquette per wall per slice ⇒ `4T` temporal plaquettes,
cross-section `=∂p₀` (the 4-link boundary of the source plaquette). Everything else costs `≥4T+1` or is
rate-separated (`≥6T`). **[T]** (`P3-FOURWALL §1.3–§1.5`)

### 1.2 Its weight, its overlap, and the modulus-square `[T]`

The bare tube weight is
```
  w(B_T) = a₀(β)\,u^{4T} = a₀\,e^{−κ t},   κ := 4\ln(1/u)/a = m₁ (strong coupling, = Gupta M_{0⁺⁺}),
  a₀(β) = (2N)^{−p}\,s ,   s ∈ ℤ_{≥1} ,  p∈{0,2} ,   a₀ ≥ (2N)^{−2} > 0 ,
```
where `s` is the product of **link recoupling multiplicities** `∫dU∏χ_{r_i}=dim Hom(1,⊗r_i)∈ℤ_{≥0}`
(Fact 0.1, `P3-DRESSED-MASS`), the source/sink `f⊗f̄→1` singlet recoupling contributing `1`. This is
the precise content of *"the atom amplitude is the modulus-square of a positive `f⊗f̄→1` overlap"*:
the amplitude `a₀` is a product of **non-negative integers and `(2N)^{-1}` factors**, so no sign or
recoupling can undo it. The **overlap** is therefore
```
  f₁ = ⟨G|O|Ω⟩ = ⟨G|φ⟩ ,     f₁² = a₀(1 − O(u))   at strong coupling   [T] (mod KL)
```
(`P3-DRESSED-MASS §3–§4`, `P3-CYLINDER-RESIDUE` for the KL/OS spectral identification `[KL]`).

> **What "the overlap" means as an operator.** The tube is *one particle* because strong-coupling
> **confinement fuses** the four source strands into a single flux tube below the two-tube threshold.
> The overlap `f₁=⟨G|φ⟩` measures how much of the local source `φ=QO|Ω⟩` fuses into that one particle.
> It is **positive and `T`-independent**; it is **not** a statistical overlap that shrinks with the
> number of tube sites.

### 1.3 The normalised overlap — the honest, scheme-free floor `[T]`

Define `M₀:=C(0)=‖φ‖²=Var(O)` and
```
  w₁ := f₁²/M₀ ∈ [0,1] ,   atom fraction = fraction of the source's spectral weight in |G⟩ .
```
Two structural facts make `w₁` the right object:
- **Dimensionless and posivity-bounded:** `0≤w₁≤1` (it is `⟨G|P_G|G⟩/‖φ‖²` of the normalised source,
  with `w₁=1` iff the source *is* the tube, `w₁=0` iff it has no weight on `|G⟩`).
- **`Z_O`-invariant `[T]`:** under the exact transport `O→Z_O O` (`G2-RESIDUE-ATOM`, `INVENT-TRANSPORT-
  INEQ`), `f₁²→Z_O²f₁²` and `M₀→Z_O²M₀`, so `w₁` is **unchanged**. R6-c showed the dimensionful floor
  `f₁²≥r₀` is scheme-dependent (only its sign is physical) and named `f₁²/M₀` as the honest ratio; here
  `w₁` *is* that ratio, and it is exactly the atom fraction. **[T]**

The **floor to prove** is `w₁≥w₀>0`, uniform in `(L,a)`, through `β=2.5`. Equivalently
`Σ_{n≥2}w_n ≤ 1−w₀` — a bound on the **leakage** of source weight into the second tube / multi-tube cut.

---

## 2 · What RP, PF and the transfer matrix can and cannot give

The KL/OS representation (`[KL]`, `P3-DRESSED-MASS §5`) gives `C(t)=∫_{[0,∞)}e^{−Et}dμ_O(E)`,
`dμ_O≥0`, with `dμ_O=Σ_n w_n δ_{E_n}+(\text{cut})`, `w_n=|⟨n|O|Ω⟩|²`, and `w₁=dμ_O({E₁})/M₀`. The three
named tools act on this object as follows.

### 2.1 Reflection positivity `[NG]`

RP ⇒ `dμ_O≥0` ⇒ `C(t)` is a **Stieltjes/completely-monotone** function: `log C(t)` convex,
`m_eff(t):=−(d/dt)log C(t)` **decreasing from above** to `m₁` (`G2D-ISOLATION-I2`, `BS-BLIND`,
`STATE-DENSITY-I2`). But these are **upper-bound-shaped** statements: a positive measure may carry
`w₁=0` (a pure two-tube **cut**) and satisfy *every* RP consequence. The free massless field is exactly
this (`w₁=0`), so any RP route to a floor would wrongly kill the free field — it is **acid-test-blind**
and returns nothing. **[NG]** for the class of RP/moment/complete-monotonicity functionals. (Witness
(C): `log C̄` is convex in *both* the `w₁>0` and `w₁=0` cases.)

### 2.2 Perron–Frobenius `[NG]`

The gauge-invariant transfer matrix `T=e^{−aH}` is positivity-improving (`PROOF-GAUSS-INVARIANT-PF`):
`T` maps the cone of non-negative invariant functions into its interior (electric heat kernel `K_t>0`
times `e^{−a\cdot magnetic}>0`). Jentzsch/Krein–Rutman gives a simple strictly-positive vacuum `Ω` and
a finite-`L` gap `λ₁<λ₀=1`. **But** the connected source `φ=QO|Ω⟩=O|Ω⟩−⟨O⟩|Ω⟩` is **sign-indefinite**
(only the *diagonal* `⟨O⟩≥0` is sign-definite). PF bounds the diagonal `⟨Ω|O|Ω⟩`, the Perron amplitude;
it says nothing about `|G⟩∈Ω^⊥`, the one face where the single RP cone degenerates. Quantitatively, the
Hilbert-metric contraction `T^nψ/‖T^nψ‖→|Ω⟩` holds for `ψ` in the cone — *not* for `φ` — and the
overlap `⟨G|φ⟩` is the difference of two positive-vector flows, which PF cannot separate from zero.
This is `INEQ-R1`'s obstruction, relocated onto `ker Gauss`; the free Maxwell field is PF with `w₁=0`.
**[NG]** (PF removes the vacuum; the edge is what survives, and PF leaves it free).

### 2.3 The transfer matrix: a **decreasing upper bound**, not a lower bound `[T]`

This is the sharpest positive statement available, and it answers the task's tube-length question.
From the spectral representation, for any `T=t/a>0`,
```
  ŵ₁(T) := e^{E₁ t}\bar C(t) = w₁ + Σ_{n≥2} w_n e^{−(E_n−E₁)t} + e^{E₁ t}\!\int_{(E₁,∞)}e^{−Et}dμ_{\rm cut}
```
with `w_n≥0`, `E_n≥E₂>E₁`. Hence:

> **Theorem T (monotone upper bound, non-certifiability). `[T]`**
> (a) `ŵ₁(T)` is **decreasing** in `T` and `ŵ₁(T) ≥ w₁` for all `T`;
> (b) `w₁ = lim_{T→∞} ŵ₁(T)`;
> (c) so **every finite-`T` overlap estimate is an upper bound**, and it converges *downward* to `w₁`.

Witness (B): `ŵ₁(T)` starts at `1.0` (`T=0`) and decreases monotonically to `w₁=0.300`. Consequently:
- **No finite tube length certifies `w₁≥w₀`** — any finite-`T` number is an over-estimate; a small
  measured value does not *prove* a small truth, but a *plateau* at `>0` is consistent with it.
- **Combined with `R3` (Stone–Weierstrass), no finite `L` certifies it either** — a thin cut abutting
  `E₁` reproduces any finite moment list.
- The **collapse** the task asks about is therefore **in the estimator `ŵ₁(T)`** (it shrinks toward `w₁`),
  **not in the object `w₁`** (which is `T`-independent). This is the exact sense in which the tube-length
  dependence is *in the wrong direction* to help: longer tubes give *smaller* estimates, never larger.

### 2.4 The one route that *would* work — and it is circular `[O]`

A lower bound on an atom given *total data* is the **residue–level-gap lemma** L2 (`GAP-G2-edge`),
written for `w₁`:
```
  w₁ = f₁²/M₀ ≥ (Λ − Ē)/(Λ − m₁²) ,    Ē := M₁/M₀ ,  M₁ = ∫ s\,dμ_O(s) ,   nonvacuous iff Λ > Ē .
```
It **consumes the second level `Λ=E₂`** (the isolation `I₂`). Witness (E): with `Λ=E₂` the bound holds
(`w₁≥0.275≤0.400`); without a second level (`Λ→m₁²⁺`) it collapses to `≤0`; with an *unsound* `Λ` above
the true `E₂` it **over-claims** (`0.638 > 0.400`). So the only exact lower bound on the overlap is
**exactly the isolation**, which is the gap one rung up (`E2-REGRESS-CIRCULARITY` L2/L3). There is no
free exact overlap lower bound. **[O]**

> **Net of §2.** The task's requested tools — RP, PF, the transfer matrix — give, respectively: RP a
> **non-floor** (log-convexity, acid-test-blind), PF a **vacuum/gap** (edge-blind), and the transfer
> matrix a **decreasing upper bound** on the overlap (no floor). The overlap **cannot be lower-bounded
> by positivity**; it needs the polymer/cluster expansion (§3) or a certified near-tube trial (§2.4,
> circular).

---

## 3 · The only proven floor: strong coupling, and its uniformities

### 3.1 The strong-coupling overlap floor `[T] (mod [KL])`

Assemble `P3-FOURWALL` (geometry) + `P3-FULLREP-KP` (full-rep convergence) + `P3-DRESSED-MASS`
(skeleton resummation) + `P3-CYLINDER-RESIDUE` (KL/OS):
```
  f₁²_{,L} = dμ_O^{(L)}({m₁^{dressed}}) = A^{dressed}(L) ≥ ¼(2N)^{−2} > 0 ,  uniformly in L (fixed a),
  m₁^{dressed} = (4\ln(1/u) − C_{dec}u + O(u²))/a ,   A^{dressed}=a₀D_fin ,  D_fin=1±O(u) ,
```
for `β<β₁=min(1/48, c/N)`. The **dressed amplitude** `A^{dressed}` is exactly the overlap
`f₁²` (`=a₀D_fin`), and it is **positive** because `a₀>0` (recoupling) and `D_fin>0` (absolute KP
convergence, *not* termwise nonnegativity — so no blanket `c_r≥0` is used). Equivalently the
normalised `w₁=A^{dressed}/M₀ ≥ w₀(β)>0`, with `M₀=C(0)` also computable in the same expansion. **[T]**
mod `[KL]`.

**Acid test passes `[T]`.** For free massless Maxwell, `a₀=0` (the source does not fuse to a singlet
tube — no confinement), so the construction returns `f₁²=w₁=0`; in 2D `tr F²` is a conserved density
and `w₁=0` exactly (`R5 §1.1`); in both the *recoupling* may be present but the *atom* is absent (§4).
The floor is non-vacuous precisely because YM recouples `f⊗f̄→1`.

### 3.2 Volume uniformity `[T]`

`P3-DRESSED-MASS §2.4`: on `T³_L` the decoration gas may include polymers wrapping a spatial cycle,
which cost `≥L/a` and weight `≤u^{L/a}=e^{−cL}`. Since the KP criterion `uK<1` is **per-plaquette** and
`L`-independent, `A^{dressed}(L)=A^{dressed}(∞)+O(e^{−cL/a})`. So the floor is **uniform in `L`**:
`inf_L A^{dressed}(L) ≥ ¼(2N)^{-2}>0`. The transverse-position sum of the backbone (`R5`'s "geometric
multiplicity") is `O(1)` and folded into `a₀`, not an `L³` entropy (it cancels in the normalised
residue). **[T]**

### 3.3 Tube-length uniformity — the decisive point `[T]/[C]`

The overlap `f₁²` is the **residue of a pole**, hence a property of the `t→∞` transfer spectrum, hence
**`T`-independent**. Concretely, in the skeleton expansion the only `T`-dependence sits in
`Z_dec(B_T)=e^{f_{dec}t}D_fin(T)` with `D_fin(T)=1+O(u)` (`P3-DRESSED-MASS §3,§4`): the decoration gas
renormalises the **mass** (`f_{dec}`) and leaves the **amplitude** `a₀D_fin` finite and `T`-uniform.
Physically, transverse **roughening** (the flux tube's wandering) is a *decoration* — a mass (rate)
shift, not a residue suppression — and in 4D it is **logarithmic** (Lüscher term `∼log T` in the
potential) `[std]/[C]`, so even in the string picture the residue softens at most as `1/(c+log T)`, not
a power-law collapse.

> **Answer to "does the bound collapse with the tube length?" `[T]`** *No.* The overlap `f₁²` (=`a₀D_fin`)
> and the normalised `w₁` are **tube-length independent**; the strong-coupling floor is `T`-uniform. What
> *does* change with `T` is the **finite-`T` estimator** `ŵ₁(T)` (§2.3), which decreases *toward* `w₁` —
> so longer tubes give *tighter (smaller) upper bounds*, which is the **wrong direction** for
> certifying a floor but is **not** a collapse of the object. The genuine tube-length danger is not
> "collapse" but **dissolution** at fixed `T` (§5): if the binding margin →0, `w₁→0` at *any* `T`.

---

## 4 · Tests: 2D, 3D, free field — isolating what the overlap does and does not detect

The overlap `w₁` is the modulus-square of a positive recoupling. The tests below separate **"has the
recoupling"** from **"has the atom"**, and locate the discriminant precisely.

### 4.1 2D pure YM — recoupling present, atom exactly zero `[T]`

In `1+1`D, `tr F²=tr E²` **is** the energy density: conserved, diagonal in the character basis, so
`⟨0|tr F²|G(R)⟩=C₂(R)⟨χ₀|χ_R⟩=0` for `R≠0` (`R5 §1.1`). The four-wall *recoupling* `f⊗f̄→1` is a
perfectly positive integer — yet `f₁²=w₁=0` **exactly**, because the connected spectral measure has no
atom at `m>0` (it is a contact in `t`). **This is the sharpest acid test specific to the overlap
claim:** a positive recoupling `f⊗f̄→1` is **necessary but not sufficient** for the atom. The overlap
*amplitude* (`f⊗f̄→1`) is not the discriminant; the discriminant is whether that fusion produces an
**isolated particle** (dynamical persistence / isolation), not just a group-theoretic channel.

### 4.2 3D pure SU(N) YM — same four-wall mechanism, atom exists `[T-standard]/[V]`

The 4-wall count is `|P|≥4T` in 3D too (a spatial plaquette has 4 boundary links regardless of spatial
dimension, and each needs a temporal partner). So the strong-coupling residue theorem P3 transcribes
verbatim: `f₁²≥r₀(β)>0` at `β<1/48`, `L`-uniform. Lattice 3D YM is confining, gapped, with a `0⁺⁺`
atom `[V]`. 3D is **easier** for the isolation conjunct (one fewer decay channel; the no-dissolution
margin `m_{0⁺⁺}<2m_{\rm lightest}` is easier, `R6 §2.5`). So 3D supports the overlap being positive
and does not add a mechanism beyond 4D. **[T-standard]/[V]**

### 4.3 Free field — `w₁=0` exactly (acid test) `[T]`

A composite `tr F²` of a free field is a **two-particle bilinear**: its KL measure is the
two-particle cut on `[4m²,∞)` (massive) or `[0,∞)` (massless), absolutely continuous with **no delta**.
So `w₁=dμ({m²})=0` for every `m>0` (`R5 §1.6`, witness (D)). The "tube" does not exist
(`a₀=0`, no confinement recoupling). **Acid test passes discriminatingly:** every tool above returns
`w₁=0` on the free field (no manufactured floor), and `>0`-capable only through the YM recoupling.

> **Reading of §4.** The overlap `w₁` is *created by* binding (the `f⊗f̄→1` fusion into an isolated
> tube) and *absent* in 2D (conserved, no isolation) and the free field (no confinement). 3D behaves like
> 4D. The 2D case is the decisive discriminator: **positive recoupling with zero atom** shows that no
> *static* "overlap" object — only the *dynamical* isolation — can carry the floor.

---

## 5 · `β=2.5`: does anything survive?

### 5.1 The strong-coupling floor is outside its radius `[T]`

The P3 floor is a Kotecký–Preiss cluster expansion; its hypothesis is `u·K<1`, `K≤CN` with `C∼6×10²`
(`P3-FULLREP-KP`: `K:=2λ_aC_ge≤80e^{2+β}N≈6.1×10²N`). Using the measured `u=I₂/I₁` (witness (F)):
```
  β        u       uK (N=2)     inside uK<1?
  1/48   1.04e-2     12.7         no (crude K; rigorous radius is β<1/48)
  1.0    4.33e-1    528           no
  2.5    7.19e-1    877           no
```
At `β=2.5`, `uK∼10²–10³≫1`: the expansion is **far outside** its convergence radius. The floor
`f₁²≥¼(2N)^{-2}` is therefore **not defined** at `β=2.5` — not merely untransported. This is the
`1/48` Ursell wall of `R5`/`PROOF-P3 §5`, quantified: the wall is not a technicality but a factor
`∼10²` in the KP parameter. **[T]** *(The crude `K` over-estimates the radius by an `O(1)` factor; the
message `uK≫1` at `2.5` is robust — even `K=O(10)` gives `uK≫1` because `u(2.5)=0.719`.)*

### 5.2 What is *not* killed: existence and physical scaling `[V]/[C]`

- **Existence `f₁²>0`, `w₁>0` is protected** by the confined single phase + `0⁺⁺`-lightest ordering +
  the trace-anomaly carrier overlap (`R6 §5`, (H-a,b,c)) — conditional, not `[T]`.
- **Physical scaling does not collapse:** dimensionally `f₁²=|⟨0|tr F²|G⟩|²∼Λ⁶`, and the `β=2.2` lattice
  measurement gives `f₁²_bare=7.14±0.42` (`17σ`, `L=6,N=4000`; `R5 §5`,`R4`) — an `O(1)` (lattice-unit)
  overlap, *not* exponential. The `β=2.5` 40-config pilot prefers a positive `C₁/C₀∼0.18–0.29` with no
  plateau (`R5 §5`). So the **overlap does not physically collapse**; it is a robust `O(Λ³)` decay
  constant. **[V]/[C]**
- **But none of this is a proven uniform lower bound**, and it does not address isolation.

### 5.3 The obstruction is not the overlap magnitude — it is isolation `[T]/[O]`

By §2.4, any *proven* lower bound on `w₁` given total data is the residue–level-gap lemma, which is a
lower bound on `Λ=E₂` = isolation. So at `β=2.5` the overlap's *magnitude* is not the hard object (it is
robust in models and in scaling); the hard object is whether the fused tube is an **isolated** discrete
state, `Λ>m₁²+δ`, rather than a threshold resonance that drains `w₁` into the two-tube cut. This is
exactly R6's split: **residue robust, isolation = gap.** T2 sharpens the *residue* to `w₁` (scheme-free)
and confirms it is neither the collapse-prone object nor the free one.

---

## 6 · Witness — `research/v3/w7/t2_tube_overlap.py`

`py -3.12 -I t2_tube_overlap.py` (exit 0; numpy 2.4.6; deterministic, ~0.2 s). Quoted output (abridged):

```
(A) w1 = weight_lightest / total = 0.300000 ,  in [0,1] ; sum_n w_n = M0
(B) w1hat(t)=e^{E1 t}Cbar(t): t=0 -> 1.000000 ; t=8 -> 0.300097 ; t>=16 -> 0.300000
    monotone non-increasing in t : True ;  limit == w1_true : True
    => finite tube length OVER-estimates the overlap (upper bound only)
(C) RP/log-convexity: w1>0 (atom+cut) convex, m_eff 1.892->1.000 ;
    w1=0 (pure cut) convex, m_eff 2.342->2.002  => RP holds in BOTH, blind to w1
(D) dissolution: delta=1e0 -> w1=0.347 ; 1e-2 -> 0.050 ; 1e-4 -> 0.005 ; 1e-6 -> 5e-4
    => w1 -> 0 continuously as binding -> 0 (measure stays RP); free field w1=0 exactly
(E) residue-level-gap: true w1=0.400 ; Lambda=E2 -> bound 0.275 (valid, loose);
    Lambda=m1^2+ -> <=0 (vacuous) ; Lambda=3.0>E2 -> 0.638 (OVERCLAIMS)
    => the only true lower bound CONSUMES Lambda=E2
(F) a0 >= (2N)^-2 = 0.0625 > 0 ; uK: 12.7 (1/48), 528 (1.0), 877 (2.5)  => 2.5 OUTSIDE radius
```
Readings: (A) `w₁` is a clean `[0,1]` projection fraction; (B) the finite-`T` estimator is a decreasing
upper bound — the tube-length "collapse" is in the estimator, not the object; (C) RP/log-convexity is
blind to `w₁` (holds for `w₁=0`); (D) `w₁→0` continuously at dissolution threshold, free field `w₁=0`;
(E) the only exact lower bound consumes `E₂`; (F) the strong-coupling floor is `uK≫1` outside its radius
at `β=2.5`. **Witness, not proof; the theorems are §1–§5, the numerics illustrate the inequality
structure and the KP-radius fact.**

---

## 7 · Ledger

| item | statement | label |
|---|---|---|
| minimal polymer | straight four-wall tube `B_T`, `4T` temporal plaquettes, `\|P\|≥4T`, equality iff straight | **[T]** (`P3-FOURWALL`) |
| tube weight / amplitude | `w(B_T)=a₀u^{4T}`, `a₀=(2N)^{-p}s>0`, `s∈ℤ_{≥1}` (singlet recoupling `f⊗f̄→1`) | **[T]** |
| overlap at strong coupling | `f₁=⟨G\|O\|Ω⟩=⟨G\|φ⟩`, `f₁²=a₀(1−O(u))`, modulus-square of positive recoupling | **[T]** (mod `[KL]`) |
| normalised overlap | `w₁=f₁²/M₀ ∈[0,1]`, atom fraction, **`Z_O`-invariant** (sharpens R6-c) | **[T]** |
| RP route | `dν≥0`, log-convex `C(t)`, `m_eff↓m₁`; **no floor**, acid-test-blind (`w₁=0` legal) | **[NG]** |
| PF route | positivity-improving, vacuum+gap; edge in `Ω^⊥`, `φ` sign-indefinite ⇒ no overlap bound | **[NG]** |
| transfer-matrix estimator | `ŵ₁(T)=e^{E₁t}C̄(t)` **decreasing upper bound**, `ŵ₁(T)↓w₁` | **[T]** |
| non-certifiability | no finite `T` (Thm T) or finite `L` (R3) certifies `w₁≥w₀`; only `t,L→∞` | **[NG]** |
| residue–level-gap | `w₁≥(Λ−Ē)/(Λ−m₁²)`, nonvacuous iff `Λ>Ē`; **consumes `Λ=E₂`=isolation** | **[T]** / **[O]** |
| volume uniformity (SC) | `A^{dressed}(L)=A(∞)+O(e^{−cL/a})`, `inf_L≥¼(2N)^{-2}>0` | **[T]** (`P3-DRESSED §2.4`) |
| tube-length uniformity | `f₁²`, `w₁` **`T`-independent** (residue of pole; roughening=decoration=mass); no collapse | **[T]**/ **[C]** (4D log) |
| dissolution (fixed `T`) | `w₁→0` continuously as binding→0 (threshold residue); free field `w₁=0` | **[T]** (witness D) |
| 2D YM | recoupling positive, `f₁²=w₁=0` exactly (`trF²=H` conserved) — **recoupling necessary, not sufficient** | **[T]** |
| 3D SU(N) | same `4T` four-wall; P3 transcribes; atom exists, isolation easier | **[T-standard]/[V]** |
| free field | composite cut ⇒ `w₁=0` exactly (acid test) | **[T]** |
| `β=2.5` SC floor | `u(2.5)K∼10²–10³≫1`: outside KP radius, floor **not defined** | **[T]** |
| `β=2.5` object | existence protected (R6); scaling `f₁²∼Λ⁶`, `β=2.2` `f₁²=7.14` ⇒ no physical collapse | **[V]/[C]** |
| **uniform overlap floor** | `w₁≥w₀>0`, `(L,a)`-uniform, `β` through `2.5` — **not proven** | **[O]** |
| **missing estimate = ?** | `w₁≥w₀` ⟺ leakage bound `Σ_{n≥2}w_n≤1−w₀` ⟺ isolation `I₂` (L2) ⟺ **`H*` = the gap** | **[O]** |
| **Clay** | unbuilt | **[O] · 0%** |

---

## 8 · The exact missing estimate and its falsification test

**Missing estimate (stated once, exactly).**
> There exist `w₀>0` and `δ>0` such that, **uniformly in `(L,a)` and for all `β` through the physical
> `β=2.5`**, the normalised flux-tube overlap obeys
> ```
>   w₁(β;L,a) = f₁²/M₀ ≥ w₀ > 0 ,     equivalently  Σ_{n≥2} w_n ≤ 1 − w₀ ,
> ```
> i.e. the source's spectral weight in the isolated four-wall tube is bounded below — the source is
> "mostly one tube" — with the **isolation margin** `E₂−m₁²≥δ` (`I₂`) that makes `Λ` usable in L2.

By L2 this is **equivalent** to isolation, which (Single-Scale Covariance, `E2-REGRESS-CIRCULARITY`
L3) is the gap up to a fixed ratio — i.e. **`H*`**. There is no weaker overlap-only sub-target that is
not already `H*` or already known (`R6`): the residue *sign* is protected (P3 strong coupling), the
residue *floor* is robust in models (double well), but the *uniform, transported* floor at `β=2.5`
**is** the isolation.

**Falsification test (the decisive, cheap experiment).**
Measure, on SU(2) `β=2.5`, `L=8,12`, the zero-momentum `0⁺⁺` `tr F²` correlator (`R4` protocol: APE family
`{0,2,4,8}`, GEVP + two-state fit), and extract
```
  w₁(L,T) = f₁²(L,T) / M₀(L)  at each T=t/a and L .
```
- **Floor HOLDS** iff `w₁(L,T)` reaches a **plateau** that is **`T`-independent** (`∂_T w₁→0`) and
  **`L`-independent** (`∂_L w₁→0`) with plateau value `w₀>0`. Predicted `w₀=O(0.2)–O(0.5)` from the
  `β=2.2` pilot (`C₁/C₀∼0.18–0.29`) and the physical `f₁²_bare=7.14`.
- **Floor FALSIFIED** iff there is a sustained **`T`-drift** `∂_T w₁<0` (leakage into the two-tube cut as
  the tube is followed longer — the true residue is *below* every finite-`T` estimate, Thm T) or a
  sustained **`L`-drift** `∂_L w₁<0` (Stone–Weierstrass thin-cut leakage). Either drift, if it does not
  vanish as `T,L→∞`, drives `w₁→0` and kills the floor.
- **Acid cross-checks:** the same estimator must return `w₁=0` on a free-field / `β=0`-Haar ensemble
  (no tube) and on 2D YM (conserved `trF²`); if it returns `>0` there, the estimator is reading a
  carrier artifact, not the atom.

The test is the *same* data as R4; the new content is that it is decisive for the **scheme-free**
`w₁`, and that the *sign of the drift* (`∂_T w₁`) is a sharp, pre-registered discriminator between
"isolated tube" (plateau) and "dissolution" (drift).

---

## 9 · Sharpest line

**The minimal connected `tr F²` polymer is a single straight four-walled flux tube; its overlap with the
lightest `0⁺⁺` is, in the only scheme-free normalisation, the atom fraction `w₁=f₁²/M₀∈[0,1]` — a
modulus-square of the positive `f⊗f̄→1` recoupling, `O(1)` and positive at strong coupling, and
*invariant* under the `Z_O` transport that R6-c showed spoils the dimensionful floor. Reflection
positivity gives only log-convexity (blind to `w₁`, holds at `w₁=0`), Perron–Frobenius gives only the
vacuum and the gap (`φ` is sign-indefinite, the edge is in `Ω^⊥`), and the transfer matrix gives only a
*decreasing upper bound* `ŵ₁(T)` on the overlap — so finite `T` (or finite `L`) can never certify a
floor, and the tube-length "collapse" is in the **estimator**, not the **object** (`w₁` is
`T`-independent; roughening is a mass shift, a logarithmic one in 4D). The one proven positive overlap
bound is strong coupling (`f₁²≥¼(2N)^{-2}`, `L`-uniform, `β<1/48`), and at `β=2.5` the KP parameter is
`uK∼10²–10³≫1` — outside its radius, so the bound is not even defined there. The object does not
physically collapse (existence protected by confinement + trace anomaly; `f₁²∼Λ⁶`, `β=2.2` data
`f₁²=7.14`), but no uniform positive lower bound is proven; by the residue–level-gap lemma any such
bound *is* isolation `I₂` at `E₂`, i.e. `H*`, i.e. the gap. The exact missing estimate is a uniform
`w₁≥w₀` (equivalently `Σ_{n≥2}w_n≤1−w₀`); its falsification is a sustained `∂_T w₁<0` or `∂_L w₁<0`
drift of the `0⁺⁺` `tr F²` atom fraction. **No QED. Clay 0%.**

---

*Refs.* `R-WAVE-SYNTHESIS.md` (§7 next species), `R5-EXACT-ANCHORS.md` (§1,§2,§5),
`R6-ATOM-ADVERSARIAL.md` (§2.3,§4,§5); `unconditional/PROOF-P3-UNCONDITIONAL.md`,
`P3-FOURWALL-DERIVATION.md` (§1–§4), `P3-FULLREP-KP.md`, `P3-DRESSED-MASS.md` (§0–§5),
`P3-CYLINDER-RESIDUE.md` (`[KL]`), `PROOF-GAUSS-INVARIANT-PF.md`, `G2D-ISOLATION-I2.md`,
`GAP-G2-edge.md` (residue–level-gap), `BS-BLIND.md`, `STATE-DENSITY-I2.md`, `E2-REGRESS-CIRCULARITY.md`
(L2,L3), `G2-RESIDUE-ATOM.md`, `INVENT-TRANSPORT-INEQ.md`; `v3/SM2-RANK-ONE.md`, `H2-DEFECT-BOUND.md`.
**Standard (textbook, cited from corpus, not re-fetched):** Kotecký–Preiss CMP **103** (1986) 491;
Osterwalder–Seiler Ann. Phys. **110** (1978) 440; Källén (1952)/Lehmann (1954); Lüscher CMP **54**
(1977) 283 (roughening/Lüscher term); Krein–Rutman/Jentzsch; Reed–Simon IV §XIII; Gupta
arXiv:hep-lat/9807028 §14.2.2 (four-wall `M_{0⁺⁺}`); Münster PLB **95** (1980) 59. **Witness:**
`research/v3/w7/t2_tube_overlap.py` (deterministic, numpy only, exit 0). **No third-party downloads; no
`[C]` upgraded to `[T]`; no QED; Clay 0%.**
