# U1-HARNACK-BOUND — a PF / Doob-h-transform / Harnack / Sturm–TP₂ lower bound on the tube overlap `w₁=f₁²/M₀`: the *vacuum* h-transform is trivial, Harnack controls the eigenfunction's **shape** but never the **source alignment**, TP₂/Sturm supplies an **order** that `𝒜/𝒢` lacks *and* that does not help anyway — so the three named tools give **no** floor `w₁≥w₀`, and the exact missing estimate is the alignment inequality `(E-U1)`, which *is* isolation `I₂` one rung up

> **2026-11-24 · Wave U1 — the missing object of the T wave (a uniform `w₁≥w₀`).**
> T2 proved the atom fraction `w₁=f₁²/M₀∈[0,1]` of the minimal four-wall tube has only *upper* bounds (RP
> log-convexity; vacuum-anchored PF; decreasing transfer-matrix bound) and named the missing object: a
> uniform `w₁≥w₀>0`, equivalently the leakage bound `Σ_{n≥2}w_n≤1−w₀`. The task here is narrower and
> sharper: use the **ground-state (Doob) transform**, the **Harnack inequality**, and **Perron-eigenvector /
> Sturm–TP₂ comparison** to convert the vacuum-anchored PF into a lower bound on the *first excited* overlap,
> and decide whether the first excited eigenvector's **node structure** can give `w₁≥w₀`. We do exactly that,
> decide it, test 1D/2D/free, and state the precise missing estimate.
> **Reads.** `T-WAVE-SYNTHESIS.md`, `T2-TUBE-OVERLAP.md`, `T6-BINDING-ADVERSARIAL.md`; `R5-EXACT-ANCHORS.md`,
> `R2-ISOLATION.md`; `T3-TP2-ISOLATION.md`; `unconditional/PROOF-P3-UNCONDITIONAL.md`,
> `PROOF-GAUSS-INVARIANT-PF.md`. **Witness.** `research/v3/w7/u1_harnack_witness.py`
> (`py -3.12 -I`, numpy+scipy, deterministic, exit 0; output quoted §7). No third-party downloads.
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional · **[NG]** no-go for a named
> class · **[O]** exact missing step · **[V]** numerical witness here · **[std]** textbook tool used as such.
> **No QED. Clay 0%.**

---

## 0 · Verdict (up front; no QED, Clay 0%)

1. **The master identity — `w₁` is a source *alignment*, not a positivity datum `[T]`.** With
   `O=¼tr F²≥0`, `⟨O⟩=m_O`, `o:=O/m_O≥0`, `E_μ[o]=1`, `φ=QO|Ω⟩=m_O(o−1)`, and `T`'s eigenbasis
   `{ψ_n}` (`ψ₀=Ω`),
   ```
       w₁ = ⟨ψ₁|φ⟩²/‖φ‖² = ⟨ψ₁|o⟩²/Var_μ(o) ,     Σ_{n≥1} w_n = 1 .                          [T]
   ```
   **`w₁` is the squared direction cosine of the source `o−1` onto the first excited eigenvector.** A
   lower bound on `w₁` is therefore a statement that *the local density `¼tr F²` is not nearly orthogonal
   to the lightest `0⁺⁺` mode.* That is a statement about `(source, eigenvector)`, not about `T`'s
   positivity.

2. **The Doob h-transform buys nothing for the vacuum and is not a Markov semigroup for `ψ₁` `[T]/[NG]`.**
   The gauge-invariant transfer matrix has **ground state `Ω≡1`** (constant), so the vacuum-anchored
   h-transform `T_h f=(λh)^{-1}T(hf)` with `h=Ω`, `λ=1` is **the identity**. For `h=ψ₁`: `ψ₁` changes
   sign, and `T_{ψ₁}` is not positivity-preserving across the node (`T` is a *global* positive kernel, and
   `T(ψ₁f)` receives negative contributions from the opposite lobe for `f≥0`). So the h-transform cannot be
   globalised across the nodal set; the minimal repair is **two separate cones on the two nodal domains with
   a relative sign** — i.e. an order-2 / second invariant cone, which `𝒜/𝒢` does not carry (R2 §3.3,
   `PROOF-GAUSS §3b`).

3. **Harnack controls the eigenfunction's *shape*, never the *source alignment* `[T]/[NG]`.**
   On a nodal domain `D⊂{ψ₁>0}`, `ψ₁` satisfies a Harnack inequality `sup_D ψ₁≤κ_D inf_D ψ₁`; this gives
   `|ψ₁|≥δ₀` on `D` up to a *shape* constant. But `w₁` involves the *product* `ψ₁·(o−1)`, and `o−1` is not
   part of the Harnack data. **The one conditional bound Harnack *can* give (Lemma U1-C, §3) multiplies a
   Harnack shape `δ₀` by a source-sign-alignment hypothesis and a source-tameness ratio `a/b`; the
   alignment is the missing input, not something Harnack supplies.**

4. **Sturm/TP₂ (node theorem) does NOT give `w₁≥w₀` `[NG]` — for two independent reasons.**
   (a) *Order is absent where it matters.* The node theorem / total positivity needs a **linear order**;
   `𝒜/𝒢` is orderless (R2 §3.3, `INVENT-SECOND-CONE`, `H6 §8`), and even the symmetric 2-space (a cycle)
   has no Sturm order. (b) *Even WITH order it fails.* A Gaussian/Mehler heat kernel is **strictly TP₂** on
   the line — nodes counted, Sturm-ordered — yet a **number-parity-even** composite overlapping a
   **number-odd** first excited state has `w₁=0` *exactly* (§7·D,E). Order fixes `ψ₁`'s shape; it does not
   fix which way the source points.

5. **The decisive no-go: `w₁` is not a functional of the positivity structure `[T]/[NG]`.** PF-positivity,
   the Harnack constant, the gap `λ₀−λ₁`, and the ratio `λ₁/λ₂` are **functions of `T` alone**. But the
   source is free in `Ω^⊥`: on a fixed `T` one can choose `o−1=ε(cosθ ψ₁+sinθ ψ₂)` and get
   `w₁=cos²θ→0` (§7·B). So **no bound that uses only PF + Harnack + gap can lower-bound `w₁`.** The free
   field realises the no-go *physically* (`A→−A` parity, §5) — it is not a toy-only statement.

6. **Test battery `[V]` (must all agree).** 1D ordered path: node count `=n` (Sturm present), conditional
   Harnack bound holds; 1D free composite: the one-particle state is parity-invisible to the composite and
   the states it reaches are the two-particle cut ⇒ `w₁=0`, while the kernel is positive *and* strictly TP₂;
   2D `T²` rotor: the first excited level is 4-fold and the perfectly local source `cos(θ₁−θ₂)` has `w₁=0`
   exactly (it sits in level 2); a point source on a delocalized level has `w₁=2/(L−1)→0`. **Free field
   fails the prediction, as required.**

7. **At `β=2.5` the bound does not survive `[O]`.** The Harnack *shape* constant `κ_D` of the 4D transfer
   matrix scales like `e^{cV}` (`PROOF-GAUSS §2`, `P4 §3`), so `δ₀→0` as `V→∞`; and the *alignment* is
   unknown (and source-dependent, and false for free). The exact missing estimate is `(E-U1)`:
   ```
       uniform in (L,a), β through 2.5:   ⟨ψ₁|o−1⟩² ≥ w₀ Var_μ(o)  (equivalently w₁≥w₀)           [O]
   ```
   plus a uniform Harnack/shape constant. By the master identity `(E-U1)` **is** `w₁≥w₀`, which is
   isolation `I₂` one rung up (`R2 Lemma D`, T2 §2.4), i.e. `H*`. **The Harnack route factors the wall into
   (tractable shape) × (alignment = the wall) and does not close it.**

**One line.** *`w₁` is the source's direction cosine onto `ψ₁`; the vacuum Doob transform is the identity,
the `ψ₁` transform is not a Markov semigroup, Harnack sees the eigenfunction's shape but not the source's
alignment, and Sturm/TP₂ supplies an order that `𝒜/𝒢` lacks and that does not help even when present (free
field: strict TP₂, `w₁=0` by parity). All three tools give only a **conditional** bound
`w₁≥(δ₀ a/b)²μ(supp)` whose non-positivity input is the source-sign alignment; the exact missing estimate
is the alignment inequality `(E-U1)`, which is isolation `I₂` = `H*`. No QED. Clay 0%.*

---

## 1 · The master identity: `w₁` is a source alignment

Fix the finite-volume physical Hilbert space `H_phys=L²(𝒜/𝒢,μ)` (μ normalized, gauge-invariant), the
RP/OS transfer matrix `T=e^{−aH}` (self-adjoint, compact, positivity-improving, `0≤T≤1`), and its
orthonormal eigenbasis `Tψ_n=λ_nψ_n`, `1=λ₀>λ₁≥λ₂≥…`, with `ψ₀=Ω≡1` (constant — the invariant measure is
the PF measure, `PROOF-GAUSS §2`). The carrier is `O=¼tr F²≥0` (a positive operator density), `⟨O⟩=m_O>0`,
`o:=O/m_O≥0`, `E_μ[o]=∫o\,dμ=1`; the connected source is `φ=QO|Ω⟩=O−⟨O⟩=m_O(o−1)`, and
`M₀=‖φ‖²=m_O²Var_μ(o)`.

> **Theorem U1-1 (master identity) `[T]`.** With `w_n:=|⟨ψ_n|φ⟩|²/M₀` and `ψ₀=1`,
> ```
>    w_n = ⟨ψ_n|o⟩²/Var_μ(o)  (n≥1) ,       w₁ = ⟨ψ₁|o⟩²/Var_μ(o) ,       Σ_{n≥1} w_n = 1 .
> ```
> *Proof.* `⟨ψ_n|φ⟩=m_O(⟨ψ_n|o⟩−⟨ψ_n|1⟩E[o])=m_O⟨ψ_n|o⟩` for `n≥1` (since `⟨ψ_n|1⟩_μ=0`). Divide by
> `M₀=m_O²Var(o)`. The weights sum to `‖P₁^⊥φ‖²/‖φ‖²=1` because `φ⊥Ω`. ∎

Three immediate consequences, which frame the whole problem:

- **`w₁` is Z_O-invariant and dimensionless** — already T2 §1.3, `[T]`; the present note inherits it.
- **`w₁` is a *direction cosine*:** `w₁=cos²∠(φ,ψ₁)` in `Ω^⊥`. It depends on the *pair*
  (source, first excited vector). `0≤w₁≤1`, `w₁=1` iff `φ∥ψ₁`, `w₁=0` iff `φ⊥ψ₁`.
- **Positivity lives entirely in `T`; alignment lives in `φ`.** Every tool named in the task (PF, Doob,
  Harnack, Perron comparison, Sturm/TP₂) is a property of `T` (its kernel/cone/order). *None* refers to
  `φ`. This is the structural reason all three fail — made precise in §2–§5.

Lemma U1-1 also separates the two facts the corpus keeps apart: the *ceiling* `χ` (a harmonic mean of `T`'s
spectrum, R5) is a `T`-quantity; the *atom* `w₁` is a `(T,φ)`-quantity. Positivity can pin `χ`; it cannot
pin `w₁`.

---

## 2 · The Doob h-transform: trivial for the vacuum, non-global for `ψ₁`

The **ground-state (Doob) transform** of a positive kernel with a positive eigenfunction `h` (`Th=λh`) is
```
    (T_h f)(x) = (λ h(x))^{-1} ∫ T(x,dy)\,h(y)\,f(y) ,        T_h 1 = 1 .                 [std]
```
It is a probability kernel that removes the `h`-tilt. The task asks whether the *vacuum-anchored* PF
(`h=Ω`) can be steered to the first excited overlap.

> **Proposition U1-2 (the vacuum transform is the identity) `[T]`.**
> For the gauge-invariant transfer matrix, `Ω≡1` and `λ₀=1`, so `T_Ω=T`. The vacuum-anchored Doob transform
> is **the identity map**: there is no tilt to flatten, because the PF cone is already *centred on the
> constant*. Consequently the h-transform has no leverage on `ψ₁`.

> **Proposition U1-3 (the `ψ₁`-transform is not a Markov semigroup) `[T]/[NG]`.**
> `ψ₁` is sign-indefinite (`⟨ψ₁|1⟩=0`, `ψ₁≢0`). Put `λ₁:=Tψ₁/ψ₁` (defined off the nodal set). For `f≥0`
> supported in one nodal domain `D_+={ψ₁>0}`:
> ```
>    (T_{ψ₁}f)(x) = (λ₁ ψ₁(x))^{-1}∫_{D_+}T(x,y)ψ₁(y)f(y)dy \;+\; (λ₁ψ₁(x))^{-1}∫_{D_-}T(x,y)\underbrace{ψ₁(y)}_{<0}f(y)dy
> ```
> the second term is **negative** for `x∈D_+` (the kernel `T>0` is *global*: it does not respect the nodal
> partition), so `T_{ψ₁}` does not map the cone of `D_+`-supported non-negative functions into itself, and
> `T_{ψ₁}` is not positivity-preserving. It is therefore not a (sub-)Markov semigroup and carries no PF
> cone. The only way to make a positivity framework is to keep **two** cones, one per nodal domain, with a
> relative sign — an **order-2** (second-invariant-cone) structure. `[T]`
> On `𝒜/𝒢` no such second cone exists (R2 §3.3; `PROOF-GAUSS §3b`; `INVENT-SECOND-CONE`). `[NG]`

**Reading.** The Doob h-transform is the textbook device that turns a ground state into a Markov process;
here it (i) does nothing at the vacuum (Ω constant) and (ii) degenerates at the first excited state (ψ₁
sign-changing, global kernel). It does not "convert the vacuum PF into a first-excited PF." The conversion
*is* the second invariant cone, and that cone is the wall — not a computation one can carry out.

---

## 3 · Harnack: shape, not alignment — and the one conditional bound it gives

The **Harnack inequality** for a positive kernel `T` and a positive (sub-)eigenfunction `u>0`,
`Tu≤λu`, bounds the ratio of `u` at two points: `u(x)≤κ(x,y)u(y)`. Applied on a nodal domain where `ψ₁`
has a definite sign, it controls `ψ₁`'s **shape**.

> **Lemma U1-4 (Harnack shape) `[T]/[std]`.** On a nodal domain `D` (`ψ₁>0` or `<0`),
> ```
>    sup_D|ψ₁| ≤ κ_D inf_D|ψ₁| ,        i.e.  |ψ₁|∈[δ₀, κ_D δ₀] on D ,  δ₀:=inf_D|ψ₁|>0 .
> ```
> This bounds the *internal oscillation* of `|ψ₁|` on `D`; it gives **no** information on (i) the relative
> L²-mass of the two lobes `p_±/‖ψ₁‖²`, (ii) the absolute normalisation, or (iii) the source `o−1`.

The reason Harnack cannot reach `w₁` is now transparent: `w₁=⟨ψ₁|o−1⟩²/‖o−1‖²` involves the **product**
`ψ₁·(o−1)`, and `o−1` is not a Harnack object. The one bound that *combines* Harnack with a source property
is:

> **Lemma U1-C (Harnack × alignment ⇒ `w₁` floor; conditional) `[T]`.**
> Let `f=o−1∈Ω^⊥`. Assume
> - (A) **sign alignment**: `f·ψ₁≥0` μ-a.e. (the source has the sign of `ψ₁`'s lobe it lives on);
> - (B) **tameness**: `a≤|f|≤b` on `supp f` (`0<a≤b`, e.g. an indicator `a=b`);
> - (H) **Harnack shape**: `|ψ₁|≥δ₀` on `supp f` (Lemma U1-4);
> - (S) **support mass**: `μ(supp f)≥m₀`.
> Then
> ```
>    w₁ ≥ (δ₀ a/b)² μ(supp f)  ≥  (δ₀ a/b)² m₀  > 0 .
> ```
> *Proof.* `⟨ψ₁|f⟩=∫ψ₁ f dμ=∫|ψ₁||f|dμ ≥ δ₀‖f‖₁ ≥ δ₀\,a\,μ(supp f)`; and `‖f‖≤b√(μ(supp f))`. Hence
> `⟨ψ₁|f⟩/‖f‖ ≥ (δ₀ a/b)√(μ(supp f))`; square. ∎ (Witness §7·C: indicator sources give
> `w₁∈{0.33,0.38,0.80} ≥ δ₀²μ` — all True; the bound is loose but real.)

**The honest accounting of Lemma U1-C.** Only **δ₀** is a positivity (Harnack) quantity. The hypotheses
**(A) alignment**, **(B) tameness**, **(S) support mass** are **source** data. In particular **(A)** — that
the local `¼tr F²` fluctuation points the *same way* as `ψ₁` on its lobe — is exactly what is *not* known
in Yang–Mills: T2 §2.2's "`φ` is sign-indefinite, PF is vacuum-anchored" is here refined to "`φ`'s sign
relative to `ψ₁` is unconstrained; if it is *anti*-aligned the bound is vacuous and `w₁` can be small." The
next section shows the alignment can even be *exactly zero* (free field, 2D rotor) or made *arbitrarily
small* (fixed `T`, rotating source).

---

## 4 · Sturm / TP₂ (node theorem): order helps neither place it is needed

The **node theorem / total positivity of order 2 (TP₂)** is the order-2 machinery: a strictly TP₂ kernel
Sturm-orders its spectrum, so `ψ_n` has exactly `n` nodes and `λ_1>λ_2` is forced (T3 §1; `[std]` Karlin,
Gantmacher–Krein). The task asks whether this sign structure gives `w₁≥w₀`. Two independent obstructions:

**(4a) Order is absent on `𝒜/𝒢` `[NG]`.** TP₂/Sturm needs a **linear order** on the base space. `𝒜/𝒢`
(the gauge quotient) carries no order (R2 §3.3, `PROOF-GAUSS §3b`, `INVENT-SECOND-CONE`); even a symmetric
cycle has no Sturm order (a cycle has no first/last point). On such a space "the first excited vector has
one node" is not a theorem, and the nodal partition is not a Sturm partition. [NG]

**(4b) Even WITH order, the node theorem does not control `w₁` `[T]/[NG]`.** The node theorem fixes `ψ₁`'s
*shape* (one sign change); combined with Harnack it fixes `|ψ₁|`'s internal comparability (Lemma U1-4). But
`w₁` is the source's projection, and a sign-changing `ψ₁` can be *perfectly orthogonal* to a fixed source for
symmetry reasons the order cannot see. Two exact witnesses:

| model | order? | TP₂? | `w₁` of the composite | why |
|---|---|---|---|---|
| 1D Gaussian/Mehler heat kernel (`e^{−tH}`, `H` oscillator) | **yes** | **strictly TP₂** (min 2×2 minor `≥0`, §7·E) | **`0`** | composite is number-*even*, first excited state number-*odd* (`A→−A`, §7·D) |
| 2D `T²` rotor (`H=−Δ/2I`, ρ²=2, exact) | — | Mehler TP₂ | **`0`** for `cos(θ₁−θ₂)` | that source sits in level `|k|²=2`; the 4-fold first level is `|k|²=1` (§7·H) |
| 1D ordered **path** (graph-Laplacian `e^{−sH}`) | yes | TP₂/Sturm (nodes `=n`, §7·G) | `∈(0,1]` **iff aligned** | with an aligned indicator source Lemma U1-C holds; with a misaligned source it is `0` (§7·C) |

**Reading.** In the one place order *is* present and one can check (the path, the Gaussian line), the node
theorem gives a *shape* — and `w₁` is still zero whenever the source is misaligned with `ψ₁`. So **the node
theorem is neither sufficient (misalignment kills it) nor statable where it matters (orderless `𝒜/𝒢`)**.

---

## 5 · The three no-go witnesses (why no positivity bound exists)

The three ways `w₁` fails to be controlled by positivity are logically distinct; all three are realised
exactly in the witness.

> **No-go U1-NG1 (source direction) `[T]+[V]`.** Fix any `T` with `dim Ω^⊥≥2` and simple `λ₁<λ₂`. For
> `φ_θ=ε(cosθ ψ₁+sinθ ψ₂)` (a legal `o=1+φ_θ≥0` for small `ε`), `w₁=cos²θ→0` as `θ→π/2`, while
> `PF`, `κ_D` (Harnack), the gap `λ₀−λ₁` and `λ₁/λ₂` are **fixed** (functions of `T` alone). Hence no
> functional of the positivity structure lower-bounds `w₁`. [V §7·B: `w₁: 1.000→0.000` at fixed
> `gap=0.2412`, `λ₁/λ₂=1.933`, `ψ₁` shape `3.73`.]

> **No-go U1-NG2 (free-field parity) `[T]+[V]`.** For free Maxwell / free scalar, `A→−A` is a symmetry
> (`F→−F`, `tr F²` invariant); the transfer matrix and vacuum are even, the **one-particle states are odd**,
> and the composite `tr F²` / `:φ²:` is **even** (number-even bilinear). So the composite **cannot reach the
> one-particle atom**: `⟨1-part|O|Ω⟩=0` *exactly* (§7·D). The states it *can* reach (number-even) start at
> the **two-particle** level, whose infinite-volume spectrum is the **pure cut** `[4m²,∞)` of R5 §1.6, with
> no delta. Hence the KL measure of the free composite is a.c. and `w₁=dμ({m₁})=0`, while `T` is
> positivity-improving, gapped (`gap=ω>0`), and even **strictly TP₂**. [V §7·D,E: `<1|:φ²:|0>=0.0e0`,
> transfer `>0`, TP₂ True; the literal "no atom" is the R5 §1.6 cut, the parity is its mechanism.] This is
> the acid test: any Harnack/PF/TP₂ theorem "proving" `w₁>0` would be refuted here.

> **No-go U1-NG3 (delocalisation) `[V]/[C]`.** A *point* source (the single-plaquette `¼tr F²` is a point
> in the dual lattice) against a **delocalised** lowest mode has `w₁~|ψ₁(x)|²/Var`, and a plane-wave mode
> has `‖ψ₁‖_∞²~1/V'.` On the `L`-cycle, `w₁=2/(L−1)→0` (§7·F, exact). A **localised** light mode (finite
> correlation length `ξ`, i.e. a genuine gap) gives `‖ψ₁‖_∞²~1/ξ^d`, whence `w₁` is *not* volume-suppressed.
> So a uniform floor needs **localisation = a mass gap with finite `ξ`**. [C: this is the physical face of
> "the lightest `0⁺⁺` is a bound state, not a plane wave".]

**Physical synthesis.** NG2 is the free-field failure *and* the sharpest statement: the obstruction is a
**representation/parity** fact (which operator sees which state), not a positivity fact. This is the U1
restatement of R5's "the discriminant is binding, not Haar/compactness/carrier": only when confinement
**binds** a `0⁺⁺` particle into the *same* symmetry channel as `O` does `⟨ψ₁|O|Ω⟩` become nonzero.
Positivity (PF/Harnack/TP₂) is orthogonal to that.

---

## 6 · Does the bound survive at `β=2.5`? — the exact missing estimate

Take stock of Lemma U1-C's four inputs at the physical point `β=2.5`, `(L,a)→` cont.:

| input | status at `β=2.5` | label |
|---|---|---|
| **(H)** Harnack shape `δ₀` on a lobe | Harnack *holds* at finite `V` (`T>0`), but the constant `κ_D~e^{cV}` (`PROOF-GAUSS §2`, `P4 §3`), so `δ₀` **not** `V`-uniform | **[O]** |
| **(A)** source sign-alignment `f·ψ₁≥0` | **unknown**; source-dependent; *false* for free (NG2); free to be made `≈0` (NG1) | **[O] = the wall** |
| **(B)** tameness `a/b` | OK for `¼tr F²` (a bounded local density) | **[T]/[V]** |
| **(S)** support mass `μ(supp f)` | OK, but point source + delocalised mode ⇒ NG3 unless localised | **[V]** |

So **the Harnack/PF/Sturm route does not produce a uniform `w₁≥w₀` at `β=2.5`.** What it *does* is factor
the wall cleanly:

> **Factorisation `[T]`.** A uniform floor reduces to
> `(uniform Harnack shape δ₀) × (source sign-alignment (A)) × (source tameness & localisation)`.
> The first factor is positivity (tractable-in-principle, blocked only by the `e^{cV}` PF constant); the
> **second is the wall** and is *exactly* the order-2 datum `𝒜/𝒢` lacks; the third is localisation =
> the gap.

> **Exact missing estimate `(E-U1)` `[O]`.** There exist `w₀>0`, `c₁>0`, `c₂>0` such that, uniformly in
> `(L,a)` and for all `β` through `2.5`, the source alignment and shape obey
> ```
>    (A′)  ⟨ψ₁|o−1⟩² ≥ c₁ Var_μ(o)        (equivalently  w₁ ≥ w₀ ,  i.e.  Σ_{n≥2} w_n ≤ 1−w₀) ,
>    (H′)  δ₀ ≥ c₂  and  ξ ≤ ξ₀  (uniform Harnack/shape + localisation) .
> ```
> By Theorem U1-1, `(A′)` **is** `w₁≥w₀`, which by R2 Lemma D is implied by isolation `I₂`
> (`Λ★−m₁²≥δ`) with weight-dominance `Λ★>R̄`; neither `(A′)` nor `(H′)` is a consequence of PF/Harnack/TP₂
> (NG1–NG3). **This is the same `H*` one rung up** — consistent with T2 §2.4/§8 and R2 §6.1.

**What would actually prove `(A′)`.** Not positivity. Two honest routes, both already named by the corpus:
- **direct dynamics** — the strong-coupling polymer *does* give alignment at `β<1/48` (P3: the minimal
  `tr F²` polymer fuses into the four-wall tube, `w₁≥w₀(β)`, `[T]` there) — but it is outside the KP radius
  at `β=2.5` (`uK≫1`, T2 §5.1);
- **the second invariant cone / binding-induced TP₂** placing `φ` in the `ψ₁`-cone with the right sign —
  the order-2 datum `(E-I2)` of R2/T3, unbuilt on `𝒜/𝒢`.

So the Harnack programme *ends exactly at the corpus's wall*, but with a sharper diagnosis: it is the
**alignment**, not the positivity, that is missing, and alignment = binding into the same channel.

---

## 7 · Witness — `research/v3/w7/u1_harnack_witness.py` (`py -3.12 -I`, exit 0)

Deterministic, numpy 2.4.6 + scipy 1.18.1, no MC, ~0.5 s. Quoted output (abridged):

```
[A] MASTER IDENTITY  w1 = <psi1|o-1>^2 / ||o-1||^2 ,  sum_n>=1 w_n = 1
    Perron vector constant-sign : True ;  gap -log(lam1/lam0) = 0.1386
    w1 = 0.062752   in [0,1] : True ;  sum_n>=1 w_n = 1.000000000000
    identity |w1 - <psi1|o-1>^2/||o-1||^2| = 0.00e+00

[B] NO-GO U1-NG1: w1 is the SOURCE direction cosine; positivity+gap do NOT bound it
    FIXED K: gap=0.2412  lam1/lam2=1.9326  psi1 shape sup/inf=3.732
     theta    w1=cos^2 th   min(o_theta)     (o>=0 legal source)
     0.000      1.000000      0.70000
     0.800      0.485400      0.70000
     1.200      0.131303      0.70000
     1.570      0.000001      0.70000     => w1 -> 0 at fixed K

[C] Harnack(shape) x sign-alignment  =>  w1 >= delta0^2 * mu(supp)
      supp=[5,7]  w1=0.3806  >= delta0^2*mu=0.0087 : True
      supp=[6,7]  w1=0.3252  >= delta0^2*mu=0.0204 : True
      supp=[5,8]  w1=0.8041  >= delta0^2*mu=0.0116 : True
    misaligned source f=psi2 : w1=9.22e-32   (=> alignment, not Harnack, is the datum)

[D] ACID TEST: free composite number-even, one-particle number-odd => w1=0
    transfer all >0 : True   gap = omega = 0.9
    <1-particle | :phi^2: |vacuum> = 0.000e+00 ;  w1(free composite) = 0.000e+00
    <2-particle | :phi^2: |vacuum> = 1.4142 (=sqrt2; the continuum)

[E] TP2/Sturm does NOT save it: Gaussian heat kernel strictly TP2 yet parity kills w1
    strict TP2 (all 2x2 minors >=0): True  (min minor 5.55e-12)

[F] LOCAL source on DELOCALIZED lowest level: w1 -> 0 as L grows (free/gapless)
     L    gap=-log(lam1/lam0)   w1 = level-1 weight   analytic 2/(L-1)
      4      2.000000          0.666667             0.666667
      8      0.585786          0.285714             0.285714
     16      0.152241          0.133333             0.133333
     32      0.038429          0.064516             0.064516
     64      0.009631          0.031746             0.031746

[G] STURM on the ordered path: psi_n has exactly n nodes  (psi_0..3: 0,1,2,3 nodes)

[H] 2D SOLVABLE (T^2 rotor, exact rho^2=2): first excited level is 4-fold
    source cos(th1)      : w1(first 4-fold level) = 1.000000
    source cos(th1+th2)  : w1(first 4-fold level) = 0.000000
    source cos(th1-th2)  : w1(first 4-fold level) = 0.000000   (a local source, EXACTLY orthogonal)
    point source: w1 = 4/(N^2-1) = 0.003910 -> 0  (delocalized mode)
```

Readings: **(A)** the master identity holds to machine zero and `Σ w_n=1`; **(B)** `w₁` slides from 1 to 0
with the source direction at *fixed* `T` (PF, gap, Harnack, `λ₁/λ₂` all constant) — the no-go U1-NG1;
**(C)** Lemma U1-C holds for aligned indicators and collapses for a misaligned source — alignment is the
datum; **(D)** the free composite is *orthogonal to the one-particle state by number parity* (it reaches
only the number-even two-particle cut), so the atom is invisible, while the kernel is positive and gapped;
**(E)** the same Gaussian kernel is *strictly TP₂*, so order does not
prevent (D); **(F)** a point source on a delocalised mode has `w₁=2/(L−1)→0` (NG3); **(G)** the ordered path
has full Sturm node structure; **(H)** in 2D the 4-fold first level is exactly orthogonal to the local source
`cos(θ₁−θ₂)`. **Witness, not proof; the theorems are §1–§6.**

---

## 8 · Model battery

| model | transfer positivity | gap | Harnack/TP₂ | first excited | `w₁` (source) | verdict |
|---|---|---|---|---|---|---|
| 1D path, `e^{−sH}`, indicator source | PF (const-signed) | `λ₀>λ₁` | Sturm, `n` nodes | one node | `>0` **iff aligned**, Lemma U1-C | tool works **conditionally** `[T]+[V]` |
| 1D Gaussian/Mehler, `:φ²:` composite | PF, `T>0` | `ω>0` | **strict TP₂** | one-particle (odd) | **`0`** (parity: reaches only the 2-particle cut) | acid PASS; order useless `[T]+[V]` |
| 2D `T²` rotor, `cos(θ₁−θ₂)` | PF, Mehler | `1/2I` | TP₂, 4-fold level | level 1 (`|k|²=1`) | **`0`** (source in level 2) | source-dependence `[T]+[V]` |
| `L`-cycle, point source | PF | `2−2cos(2π/L)→0` | TP₂ | delocalised | `2/(L−1)→0` | NG3 `[V]` |
| fixed `K`, `φ_θ` sweep | PF | fixed | fixed | `ψ₁,ψ₂` | `cos²θ→0` | no-go U1-NG1 `[T]+[V]` |
| free massless `tr F²` / `:φ²:` | PF | RP | TP₂/Maxwell | cut (no atom) | `0` | acid PASS `[T]+[V]` |

---

## 9 · Ledger

| item | statement | label |
|---|---|---|
| master identity | `w₁=⟨ψ₁\|o⟩²/Var(o)`, `Σ_{n≥1}w_n=1` — `w₁` is a source direction cosine | **[T]** |
| vacuum Doob transform | `Ω≡1` ⇒ `T_Ω=T` (identity); no tilt to exploit | **[T]** |
| `ψ₁` Doob transform | not a Markov semigroup (global kernel ignores nodal partition) | **[T]/[NG]** |
| Harnack | controls `\|ψ₁\|` shape on a nodal domain (`[δ₀,κ_Dδ₀]`), **not** the source or lobe mass | **[T]/[std]** |
| Lemma U1-C | `w₁≥(δ₀ a/b)²μ(supp)` **iff sign-aligned** source; only `δ₀` is positivity | **[T] (conditional)** |
| Sturm/TP₂, order present | fixes `ψ₁` shape / nodes; does **not** fix source alignment | **[T]** |
| Sturm/TP₂, `𝒜/𝒢` | no linear order ⇒ no second cone on the quotient | **[NG]** |
| no-go U1-NG1 | `w₁=cos²θ→0` at fixed `T` (positivity+gap+Harnack fixed) | **[T]+[V]** |
| no-go U1-NG2 | free composite number-even vs one-particle number-odd ⇒ `w₁=0`; kernel positive & TP₂ | **[T]+[V]** |
| no-go U1-NG3 | point source on delocalised level: `w₁=2/(L−1)→0`; floor needs localisation (finite ξ) | **[V]/[C]** |
| 2D rotor | 4-fold first level; local `cos(θ₁−θ₂)` exactly orthogonal ⇒ `w₁=0` | **[T]+[V]** |
| `β=2.5` uniformity | Harnack constant `~e^{cV}`, alignment unknown ⇒ no uniform `w₁≥w₀` | **[O]** |
| **missing estimate (E-U1)** | uniform `⟨ψ₁\|o−1⟩²≥c₁Var(o)` (i.e. `w₁≥w₀`) + uniform shape `δ₀≥c₂` | **[O]** |
| equivalence | `(E-U1)` = `w₁≥w₀` = isolation `I₂` (R2 Lemma D) = leakage bound = `H*` | **[T]/[O]** |
| **Clay mass gap** | unbuilt | **[O] · 0%** |

---

## 10 · Sharpest line

**`w₁=f₁²/M₀` is not a positivity quantity — it is the squared direction cosine of the local source
`o−1` onto the first excited eigenvector `ψ₁` (`w₁=⟨ψ₁|o−1⟩²/Var(o)`, both weights sum to 1). Every tool the
task names is a property of the *transfer kernel* `T`, and none refers to the *source*: the vacuum Doob
h-transform is literally the identity (the ground state is the constant, so there is no tilt to flatten);
the `ψ₁` h-transform is not a Markov semigroup (a *global* positive kernel does not respect `ψ₁`'s nodal
partition, so one needs two cones with a relative sign — the absent second invariant cone of `𝒜/𝒢`); the
Harnack inequality controls `|ψ₁|`'s shape on a nodal domain but not the source alignment; and Sturm/TP₂
supplies a node-counting order that `𝒜/𝒢` lacks *and* that does not help even where present — a Gaussian
heat kernel is *strictly TP₂* with a clean gap, yet a number-even composite is parity-orthogonal to the
number-odd one-particle atom and reaches only the two-particle cut, so `w₁=0` (the free-field acid test).
The three failure modes are distinct and all
realised: source direction (`cos²θ→0` at fixed `T`), parity (free), delocalisation (`2/(L−1)→0`). The one
bound Harnack *can* give, `w₁≥(δ₀ a/b)²μ(supp f)`, holds **only** under a source sign-alignment hypothesis
that positivity does not supply; its non-positivity input is the missing datum. At `β=2.5` the Harnack
constant is `e^{cV}` (not `V`-uniform) and the alignment is unknown, so no uniform floor survives — the
exact missing estimate is the alignment inequality `(E-U1): ⟨ψ₁|o−1⟩²≥c₁Var(o)` uniformly, which by the
master identity *is* `w₁≥w₀`, is isolation `I₂`, and is therefore `H*` one rung up. **The Harnack route
does not close the gap; it factors it into a tractable shape times the alignment wall — and the alignment
wall is binding into the same channel. No QED. Clay 0%.**

---

*Refs.* `T-WAVE-SYNTHESIS.md` (§0 verdicts, §5 next object), `T2-TUBE-OVERLAP.md` (§1.3 `w₁`, §2.2 PF [NG],
§2.4 residue–level gap, §8 missing estimate), `T6-BINDING-ADVERSARIAL.md` (§2 lattice, §7 FB1–FB6),
`R5-EXACT-ANCHORS.md` (§1.6 free-field acid, §4 `(T★)`), `R2-ISOLATION.md` (§2 Lemma A/C/D, §3.3 order tier,
§6.1 `(E-I2)`), `T3-TP2-ISOLATION.md` (§1 order-1/order-2, §2 exact ratios); `unconditional/PROOF-P3-UNCONDITIONAL.md`
(strong-coupling residue, four-wall tube, `β<1/48`), `PROOF-GAUSS-INVARIANT-PF.md` (§2 PF constants `~e^{−cV}`,
§3b edge obstruction). **Standard tools as such:** Doob h-transform (Doob 1957; Revuz–Yor *Continuous
Martingales* XII); Harnack inequality for positive kernels / eigenfunctions (Davies, *Heat Kernels and
Spectral Theory*, §3; Moser 1964; Saloff-Coste); Perron–Frobenius/Krein–Rutman (Jentzsch; Schaefer);
total positivity & Sturm oscillation (Karlin, *Total Positivity*; Gantmacher–Krein; Schoenberg–Whitney).
**Witness:** `research/v3/w7/u1_harnack_witness.py` (`py -3.12 -I`, numpy 2.4.6/scipy 1.18.1, deterministic,
exit 0, ~0.5 s). **No third-party downloads; no `[C]` upgraded to `[T]`; no QED; Clay 0%.***
