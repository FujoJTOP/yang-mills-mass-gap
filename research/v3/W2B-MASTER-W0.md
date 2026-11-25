# W2B-MASTER-W0 — the master field `W₀` **is built** in the two truncations where the question closes: the **Gaussian/Hartree** master field (exact) and the **2D Singer** field (exact), and in **both** the `tr F²` atom fraction is **exactly 0** (`g₁/Σg=0`, at every `β` incl. 2.5 — the free-field acid test at finite coupling); the **rigid-string** `W₀` has genuine binding (`m₁=2π√(2σκ)`) and gives a **positive but source-dependent** `g₁/Σg` (an **upper** estimate, not a certified floor); hence at `β=2.5` **no controlled truncation certifies a positive floor** — the floor *is* the non-Gaussian binding `= H1`

> **2026-11-24 · Wave W2B — building the vacuum master field `W₀` and evaluating `g₁/Σg` in controlled truncations.**
> W2 left `W₀` **unbuilt** in 4D and `g₁/Σg(λ=0)=0` (`W2-MASTER-FIELD.md` §4). W2B does the constructive job
> W2 deferred: it **builds `W₀` explicitly** in the three truncations the charter offers (2D exact, the
> Gaussian/Hartree large-`N` closure, and a variational large-`N` rigid-string ansatz) and **evaluates the
> atom fraction** `g₁/Σg` in each, to decide whether it has a positive floor at `β=2.5` / finite 't Hooft
> coupling. **Reads.** `W2-MASTER-FIELD.md`, `W-WAVE-CHARTER.md`, `U3-LARGE-N-BINDING.md`, `W1-STRONG-ALIGN.md`,
> `W3-SMALL-VOLUME.md`; `unconditional/CONSTRUCT-MASTERFIELD.md`, `PNT-LARGEN-LOOP.md`, `KAPPA-LARGEN.md`,
> `SLOT-DSE-FIXPOINT`/`CONSTRUCT-MASTERFIELD` §1 (Gaussian closure = DSE wall), `PROOF-EFFECTIVE-STRING`,
> `TRANSMUTE-STRING-SPECTRUM`; `R5-EXACT-ANCHORS.md` §1.1 (2D `tr F²=H`), `T2-TUBE-OVERLAP.md` §8,
> `WG-INVARIANT-ESTIMATOR.md` §5 (ceiling `χ_O`), `Q2-LOOP-BOOTSTRAP.md` §2. **Witness.**
> `research/v3/w7/w2b_master_w0.py` (`py -3.12 -I`, numpy, deterministic, **10/10 PASS**; output quoted §7).
> No third-party downloads. **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional/model ·
> **[NG]** no-go for a named class · **[O]** exact missing step · **[V]** numerical witness here · **[std]**
> textbook tool. **No QED. Clay 0%.**

---

## 0 · Verdict (up front; no QED, Clay 0%)

1. **`W₀` is built exactly in the Gaussian/Hartree truncation, and there `g₁/Σg = 0` at every `β` `[T]`/`[NG]`.**
   The large-`N` vacuum master field in the Gaussian (Hartree) closure is a **free massive functional**
   `W₀(C)=⟨exp∮A⟩_Gaussian`, i.e. a self-consistent massive propagator `G(k)=1/(k²+μ²)` with `μ=μ(β)` from its
   own tadpole (`CONSTRUCT-MASTERFIELD` §1 closure (i)). Because `O=tr F²` is **bilinear** in the gauge field,
   **Wick's theorem** makes its connected two-point function a **two-gluon bubble** `⟨OO⟩_c = 2G²`. The
   Källén–Lehmann measure is therefore **absolutely continuous** on `s≥(2μ)²` (van Hove `∼(s−4μ²)^{−1/2}`),
   and the lightest point carries **zero point mass**:
   ```
       dν(s) = a.c. ,   ν({m₁²}) = 0 ,   g₁ = ν({m₁²}) = 0 ,   g₁/Σg = 0   EXACTLY,  ∀β (incl. 2.5).   (GAUSS-0)
   ```
   This is **the free-field acid test promoted to finite coupling**: the Gaussian master field is the
   free field dressed by a tadpole, and a free (Gaussian) state has **no bound state** in a bilinear channel
   — that *is* Wick's theorem. `[T]` for the identity; `[NG]` for the class "Gaussian `W₀` ⇒ `g₁/Σg>0`".

2. **`W₀` is built exactly in 2D (Singer), and `g₁/Σg = 0` for `tr F²` but positive for a generic source `[T]`.**
   `W₀(C)=e^{−λA(C)/2}` exactly, `tr F² = H` conserved ⇒ `⟨0|tr F²|R⟩ = C₂(R)⟨χ₀|χ_R⟩ = 0` ⇒ `g₁/Σg = 0` at
   every `N` (`R5` §1.1, witness B). **But** a *generic* non-conserved 2D source (holonomy `χ_f`, or `χ_f+χ_adj`)
   gives `g₁/Σg = 1` or `1/2` — **positive and exactly computable** (witness B3). So the `tr F²` zero in 2D is a
   **selection-rule accent** (conservation), *not* a generic property of "any closed `W₀`": a solvable `W₀` *can*
   have a positive atom. This sharpens W2's 2D anchor into a positive control.

3. **The rigid-string `W₀` has genuine binding, gives `g₁/Σg>0`, but the value is source-dependent `[C]`/`[O]`.**
   The only `W₀` built so far that **binds** (a stable single particle below the two-gluon threshold) is the
   rigid surface-sum `W₀(C)=Σ_surfaces e^{−σA−κ∮K²}` (`CONSTRUCT-MASTERFIELD` closure (iv)). Its neutral `0⁺⁺`
   closed tube has `m₁=2π√(2σκ)` (`[C]`) — and its **breathing-mode variational truncation** has a discrete
   spectrum with `g₁=|⟨ψ₁|O|Ω⟩|²>0`, so `g₁/Σg>0` (witness C: `0.89–0.98` in the one-tube sector). **But the
   value moves with the identification of the source** `O=tr F²` with the tube functional (`O=V(ρ)` vs `O=ρ`
   differ by `∼0.05`), and the one-tube denominator **omits the multi-tube continuum** — so it is an **upper
   estimate of a source-ambiguous quantity**, not a certified floor. The source–tube coupling is the missing
   input `[O]` (as W2 §2.3 anticipated).

4. **Determination at `β=2.5` / finite 't Hooft coupling: no controlled truncation certifies a positive floor `[O]`.**
   The truncations split cleanly by whether they **bind**:
   ```
        NO binding   :  Gaussian W₀ (any β), 2D tr F², free field, β=0 Haar     ⇒ g₁/Σg = 0   EXACTLY
        binding      :  strong coupling (β<1/48), rigid string (any β)          ⇒ g₁/Σg > 0   but …
   ```
   The binding truncations that give `>0` either live **outside** `β=2.5` (strong coupling is inside its cluster
   radius only for `β<1/48`; `W1` §4: `uK≫1` at `2.5`) or are **source-ambiguous models** (rigid string). So the
   certified statement at `β=2.5` is a **bracket**
   ```
       0  ≤  g₁/Σg(β=2.5)  ≤  ceiling   (ceiling from the WG gauge-invariant χ_O, §5)  ,
   ```
   with the **lower face open**. The two-sided floor `g₁/Σg≥w₀>0` is **exactly** the non-Gaussian binding past
   the Gaussian truncation — i.e. `H1` in the `N`-coordinate (`W2` §6). **Building a positive `g₁/Σg` from `W₀`
   and solving the binding are the same problem.**

5. **What we do determine, quantitatively, and its honest label.** (i) The **Gaussian no-go** `(GAUSS-0)` is a
   theorem — the sharpest *controlled* `W₀` at `β=2.5` gives exactly `0`, so any claim of a Gaussian-route floor
   is acid-forbidden. (ii) The **2D positive control** (holonomy `1`, mixed `1/2`) is exact. (iii) The
   **rigid-string upper estimate** is `g₁/Σg(one-tube)≈0.9`, positive for every `κ∈[0.10,0.50]`, dropping to the
   physical `T2` band `O(0.2–0.5)` once the multi-tube leakage (`Σ_{n≥2}w_n`) is restored `[C]`. (iv) The
   **physical** floor is expected (`R5` `f₁²_bare=7.14` at `β=2.2`; trace anomaly; confinement) but **no
   truncation in hand transports it to `β=2.5`** `[O]`.

**One line.** *W2B builds the master field the W2 note deferred: the Gaussian/Hartree `W₀` is an exact free
massive functional whose bilinear `tr F²` source has a two-gluon bubble (no atom, `g₁/Σg=0` at every `β`,
Wick — the acid test at finite coupling `[T]/[NG]`); the 2D Singer `W₀` gives `0` for the conserved `tr F²` but
`1`/`1/2` for a generic source (a positive control, `[T]`); and the rigid-string `W₀` has the binding and a
positive `g₁/Σg` but its value is source-dependent (an upper estimate, `[C]/[O]`). The determination: at
`β=2.5` **no controlled truncation certifies a positive floor** — the Gaussian is exactly `0`, the binding
models are either out of cluster radius or source-ambiguous — so the floor `g₁/Σg≥w₀` **is** the non-Gaussian
binding `= H1`, and the truncations bracket it `0 ≤ g₁/Σg ≤ ceiling`. No QED. Clay 0%.*

---

## 1 · Set-up: what W2B adds, and the one object it evaluates

`W2` (`W2-MASTER-FIELD.md`) established the **ratio identity** `g₁/Σg = 1/(1+R) = ν({m₁²})/ν([0,∞))` — the
atom weight of the normalised KL measure — and that `W₀`, the equation's only input, was **unconstructed in 4D**
`[O]`. The charter (`W-WAVE-CHARTER.md` §1 W2) asks to *construct* `W₀` in a controlled truncation and compute
`g₁/Σg`. W2B builds `W₀` in the three truncations available and **evaluates the one object**

```
   g₁/Σg  =  |⟨ψ₁|O|Ω⟩|² / Var(O)  =  ν({m₁²}) / ν([0,∞))  ,      O = (1/N)tr F²  (or plaquette) ,   (§1)
```

which (W2-1) is a point mass of a positive measure. The whole W2B question is whether that point mass is
positive at `β=2.5` in a truncation that is *both* controlled (passes the acid tests) *and* at `β=2.5`.

**The three truncations, and what they are.**
- **T-I Gaussian/Hartree.** `W₀` = the self-consistent Gaussian (free massive) master field — the saddle of the
  quadratic large-`N` effective action. This is the closure `CONSTRUCT-MASTERFIELD` §1 calls (i); it is the one
  that **closes in closed form at every `β`** and reduces to the Dyson–Schwinger gap equation (`SLOT-DSE`).
- **T-II 2D exact.** `W₀ = e^{−λA/2}` (Singer), the only *exactly interacting* master field known.
- **T-III Rigid string.** `W₀ = Σ_surfaces e^{−σA−κ∮K²}` — a variational large-`N` ansatz with **binding**.

The rest of the note computes `(§1)` in each and reads off the verdict.

---

## 2 · Truncation I — the Gaussian/Hartree master field: `g₁/Σg = 0` exactly `[T]`/`[NG]`

### 2.1 The construction (closed at every `β`)

The large-`N` Gaussian master field is defined by a self-consistent two-point function. Write the gauge-field
fluctuation around the vacuum as a Gaussian with propagator `G`; the master field is
```
   W₀(C) = ⟨ e^{∮_C A} ⟩_{G} = exp( ½ ∮_C ∮_C G ),     G(k) = 1/(k² + μ²) ,     (GAUSS-W0)
```
with the mass `μ` fixed by the tadpole (the one-loop self-energy = the DSE gap equation). Two facts are exact:

- **It closes, and it exists at every `β`.** `(GAUSS-W0)` is a genuine solution of the *quadratic* truncation
  of `(MM∞)`; `μ(β)` runs from `μ(0)>0` (strong coupling) down through the weak-coupling renormalon-free mass
  (`ROOF-LARGEN-MASTER`). At `β=2.5` a massive `μ>0` Gaussian saddle exists (the DSE wall is that the *gapless*
  `μ=0` solution also persists, `CONSTRUCT-MASTERFIELD` §1/§3) — **either way `W₀` is a Gaussian**.
- **The operator `ℳ` is the free massive fluctuation operator** `ℳ = −D² + μ² (+ 2ad F_B)`, whose spectrum is
  the two-particle continuum `{2μ, 2μ + …}` — **no isolated bound state** (a free theory has none).

### 2.2 The atom fraction of a bilinear in a Gaussian state — `(GAUSS-0)` `[T]`

> **Theorem W2B-1 (Gaussian `W₀` ⇒ `g₁/Σg = 0`) `[T]`.**
> For any Gaussian master field `(GAUSS-W0)`, and any **bilinear** source `O=(1/N)tr F²` (or any operator
> bilinear in the Gaussian field), the connected source correlator is the **two-gluon bubble**
> `⟨O(x)O(0)⟩_c = 2G(x)²` (Wick). Hence its KL measure `dν` is **absolutely continuous** on `s ≥ (2μ)²` and
> `g₁ = ν({m₁²}) = 0`; therefore `g₁/Σg = 0` **for every `β`, including `β=2.5`**. ∎

**Proof.** A Gaussian measure is determined by its two-point function and satisfies Wick's theorem: all
connected correlators of order `>2` vanish, and the connected 4-point function of the *bilinear* `O∝φ²`
factorises into products of `G=⟨φφ⟩`. The Fourier transform of `2G²` is the convolution `(G∗G)(p)`, which has
a **branch cut** from the two-particle threshold `p²=4μ²` and above, and **no delta function** (a pole would need
a `1PI` four-point insertion, absent in a Gaussian measure). An absolutely continuous measure assigns zero mass
to any point, so `ν({m₁²})=0`. Witness A: the weight of the lightest spectral bin decays as `∼(nbins)^{−1/2}`
(a.c. van Hove) — `7.3e−2, 2.3e−2, 7.1e−3, 2.2e−3` at `nbins=10²–10⁵` — and the exact Gaussian Wick ratio
`Var(φ²)/Var(φ)² = 2.000` forbids a bound-state residue. ∎

**Reading.** `(GAUSS-0)` is **the free-field acid test raised to finite coupling**: the Gaussian master field is
the free field plus a tadpole mass, and a free state has no atom in a bilinear channel. It is therefore a
**no-go `[NG]` for the entire Gaussian class**: *no choice of self-consistent mass `μ(β)` — and hence no value
of `β`, `2.5` included — makes `g₁/Σg` nonzero in the Gaussian truncation.* This is the sharpest, most
controlled thing one can say with a closed `W₀` at `β=2.5`, and it is **negative for a floor**.

**Corollary (why the Gaussian truncation cannot answer the question).** The atom is **intrinsically
non-Gaussian**: it requires a connected 4-point (`1PI`) insertion in the `0⁺⁺` channel — the *binding* — which
the Gaussian closure sets to zero by construction. Resolving the floor with `W₀` **requires** going past
Gaussian, into Truncations II/III.

---

## 3 · Truncation II — the 2D exact master field: `tr F²` gives `0`, a generic source gives `>0` `[T]`

### 3.1 Singer `W₀` and the conserved-density zero

In 2D, `W₀(C)=e^{−λA(C)/2}` (Singer; `PNT-LARGEN-LOOP` §2), `λ=g²N`, `A(C)` the signed area. The spectrum is
discrete/character-labelled, but `tr F²=tr E²=H` **is the energy density**: conserved and diagonal in the
character basis, so
```
   ⟨0|tr F²|R⟩ = C₂(R)⟨χ₀|χ_R⟩ = C₂(R)·δ_{0R} = 0   for R≠0  ⇒  g₁=0  ⇒  g₁/Σg=0,   ∀N .   (2D-TRF2)
```
Witness B: SU(2) character orthogonality `max|∫χ_aχ_b−δ|=2.2e−16`, and the `tr F²` residues onto `R=0,…,4` are
`[0,0,0,0,~1e−16]` — the atom is **absent at every `N` including `N=∞`** (`W2` witness B2, `R5` §1.1).

### 3.2 The positive control — a generic source *does* give `g₁/Σg>0`

The 2D zero is a **selection-rule accent, not a generic property of a closed `W₀`**. A non-conserved source with
character content `O=Σ_R o_R χ_R` has KL weights `g_R=|o_R|²` on levels of mass `∝C₂(R)`, so
```
   g₁/Σg = |o_{R₁}|² / Σ_{R≠0}|o_R|² ,     R₁ = fundamental (lightest) ,   (2D-GEN)
```
which is **positive and exactly computable** whenever `o_{R₁}≠0`. Witness B3: holonomy `χ_f` ⇒ `g₁/Σg=1`;
`χ_f+χ_adj` ⇒ `1/2` (to `1e−12`). **This is the constructive content of Truncation II**: an exactly solvable
`W₀` *can* carry a positive atom fraction; the `tr F²` channel is special because its source is the Hamiltonian
density. It calibrates `(§1)` positively and shows the machinery works — without supplying a 4D `tr F²` floor.

> **Net of §2–§3.** Both **closed** master fields (Gaussian at any `β`; 2D exact) give `g₁/Σg=0` **for `tr F²`**,
> and the 2D model shows this is a **conservation/selection** fact, not a generic "closed `W₀`" fact. Neither
> supplies a floor at `β=2.5`; the Gaussian one forbids it within its class.

---

## 4 · Truncation III — the rigid-string `W₀`: binding, `m₁` explicit, `g₁/Σg>0` but source-dependent `[C]`/`[O]`

### 4.1 The binding and its mass

The rigid surface-sum `W₀(C)=Σ_surfaces e^{−σA−κ∮K² ds}` (`CONSTRUCT-MASTERFIELD` §1 closure (iv)) is the
**minimal variational `W₀` that binds**: a neutral (`N`-ality 0), `0⁺⁺` closed circular tube of radius `ρ` has
energy `E(ρ)=2πσρ+πκ/ρ`, with an **interior minimum** at `ρ*=√(κ/2σ)` that pure Nambu–Goto provably lacks
(`PROOF-EFFECTIVE-STRING`). The `0⁺⁺` mass is
```
   m₁ = E_min = 2π√(2σκ) ,      m₁/√σ = 2π√(2κ) ,     κ ∈ [0.138, 0.193] ⇒ m₁/√σ ∈ [3.30, 3.90] ,   (RIGID)
```
the lattice glueball band (witness C3). This is a **genuine bound state** below the two-gluon threshold — the
`[T]`-side of the `R5` discriminant ("binding, not Haar, makes the atom"). So the rigid-string `W₀` **has the
atom the Gaussian lacks**.

### 4.2 The atom fraction and why it is not a certified floor

The lightest `0⁺⁺` is the **lowest excited breathing mode** `ψ₁` of the radial Hamiltonian
`H=−(1/2M)∂_ρ²+2πσρ+πκ/ρ`. With the source `O=tr F²` identified with the tube's energy density
(`O=V(ρ)`) or with `ρ`, the KL weights are `g_n=|⟨ψ_n|O|ψ₀⟩|²`, and
```
   g₁/Σg = |⟨ψ₁|O|ψ₀⟩|² / Var₀(O)   (>0  in the one-tube / discrete sector) .   (RIGID-ATOM)
```
Witness C solves the 1D breathing Schrödinger equation (finite differences, `κ∈[0.10,0.50]`) and finds
`g₁/Σg = 0.89–0.98` (`O=V`) and `0.89–0.92` (`O=ρ`) — **strictly positive for every `κ`** (binding present). **But
it is not a certified floor**, for two independent reasons:

- **(a) Source dependence `[O]`.** `g₁/Σg` **moves with the identification** of `tr F²` with the tube functional
  (`O=V` vs `O=ρ` differ by `∼0.05`, witness C2). The `tr F²`↔tube coupling is **not fixed by the rigid-string
  functional** — exactly the missing input W2 §2.3 flagged ("gives `m₁`, not `g₁`"). The atom fraction of a
  *model* whose source is undefined is not a physical number.
- **(b) It is an upper estimate `[C]`.** The discrete one-tube denominator `Σ_{n≥1}^{tube}g_n` **omits the
  multi-tube continuum** (states at `≈2m₁`, `1/N²`-suppressed but present). The full `Σg` includes it, so the
  full `g₁/Σg ≤` the one-tube value. The physical `T2` prediction `w₀=O(0.2)–O(0.5)` is exactly the
  one-tube value minus that leakage — consistent with `0.9` shrinking to `O(0.2–0.5)`.

So Truncation III establishes **`g₁/Σg>0` structurally** (binding ⇒ the lightest mode carries weight) but
**cannot pin its value**; the value needs the source–tube coupling *and* the multi-tube leakage (isolation
`I₂`). It is a **positive upper estimate of a source-ambiguous quantity**, labelled `[C]`/`[O]`.

---

## 5 · Determination at `β=2.5`: the bracket, the ceiling, and where the floor lives

### 5.1 The bracket `[T]`

Combining the truncations and the W2 ratio identity, the certified statement at `β=2.5` is a **bracket**:
```
   0  ≤  g₁/Σg(β=2.5)  ≤  m₁²·χ_O / M₀  ≤  ceiling ,                                       (BRACKET)
```
- **lower face `= 0`**: attained *exactly* by the Gaussian `W₀` `(GAUSS-0)` and by 2D `tr F²` `(2D-TRF2)`; no
  controlled truncation pushes it above `0` at `β=2.5`.
- **upper face**: the squeeze `f₁²/χ_O ≤ m₁² ≤ M₀/χ_O` (`ch02 eq:squeeze`) gives
  `g₁/Σg = f₁²/M₀ ≤ m₁²χ_O/M₀`; with the **WG gauge-invariant ceiling** `χ_O = 0.33–0.41` at `β=2.5`
  (`WG-INVARIANT-ESTIMATOR` §5) and the physical `m₁,M₀`, this is a finite `O(1)` upper bound (the *ceiling*,
  `WG` §6). The upper face is *not* the issue.

### 5.2 Where the floor lives `[O]`/`[NG]`

The **lower face being `0`** is the whole result: it is exactly the `R5` discriminant. The atom exists iff the
dynamics **binds** `tr F²` into one particle below the two-gluon threshold — the connected `1PI` `0⁺⁺` kernel
that the Gaussian truncation sets to zero. The truncations that *do* bind are either **outside `β=2.5`** (strong
coupling `s/2>0` only for `β<1/48`; `W1` §4 `uK≫1` at `2.5`) or **source-ambiguous models** (rigid string,
§4.2). And no **moment/positivity functional** can lift the lower face, because `g₁/Σg` is a **point mass** of a
measure — identical `(M₀,M₁,M₂)` carry `g₁/Σg=0.400` vs `0.001` (witness D; `W2` §6, `[NG]`). Hence:

> **Determination (W2B).** At `β=2.5` / finite 't Hooft coupling, **no controlled truncation of `W₀` certifies a
> positive lower bound on `g₁/Σg`.** The Gaussian `W₀` — the only truncation that closes *at* `2.5` — gives
> **exactly `0`** `[T]/[NG]`. The floor `g₁/Σg≥w₀>0`, if it exists, is **exactly the non-Gaussian binding
> `= H1`**; the truncations *bracket* it (`BRACKET`) and *localise* it (binding + source-coupling + isolation),
> but do not supply it `[O]`.

### 5.3 The physical expectation (and its label)

Physically the floor is expected to be **positive**: the `β=2.2` exact anchor `f₁²_bare=7.14` (`R5` §5), the
single-dilaton trace-anomaly coupling `g₁∝m_G⁴f_G²` (`W2` §3.2, `[C]`), and confinement of the `0⁺⁺` channel all
indicate `g₁/Σg` does **not** collapse. The `T2` two-state prediction is `w₀=O(0.2)–O(0.5)` at `β=2.5`. But
this is **`[C]`/`[V-anchor]`**, not a truncation theorem: no `W₀` in hand transports the `β=2.2`/strong-coupling
value to `β=2.5` uniformly in `(L,a)`. **The floor is physically credible and truncation-open.**

---

## 6 · Acid-test battery (charter §2) — the truncations behave honestly

| # | acid test | Gaussian (T-I) | 2D exact (T-II) | rigid string (T-III) |
|---|---|---|---|---|
| 1 | **free massless** `g₁/Σg=0` | `0` exactly (Wick) `[T]` | (free limit) `0` | `0` as `κ→0`, `m₁→0` |
| 2 | **2D `tr F²`** `g₁/Σg=0` | — | `0` exactly (conserved) `[T]` | — |
| 3 | **strong coupling** `s/2` | `0` (misses the atom) `[NG]` | — | recovered in `P3` limit (out of scope) |
| 4 | **large `N`** `O(N⁰)` | `0` (exactly `N⁰`) | `1`,`1/2` exact `[T]` | `O(N⁰)` (model) |
| 5 | **continuum** `γ=0` | ratio `μ`-free `[T]` | ratio invariant `[T]` | ratio invariant `[C]` |

The **Gaussian truncation passes the free/2D zero for the right reason and *fails* the strong-coupling `s/2`** —
which is the honest signature that it truncates the binding away (it cannot have `s/2>0` because it has no atom
at all). The 2D truncation passes both `tr F²=0` and (for a generic source) a positive value. No truncation is
manufactured positive where the acid test forbids it.

---

## 7 · Witness — `research/v3/w7/w2b_master_w0.py`

`py -3.12 -I w2b_master_w0.py` (numpy 2.4.6; deterministic `seed 0`; exit 0; **10/10 PASS**; output
`w2b_master_w0.out`). Quoted output (abridged):

```
PART A  Gaussian/Hartree master field W0 : bilinear O=trF^2 has NO atom (g1/Sg=0)
   2-gluon continuum : weight in lowest bin @ nbins=1e2,1e3,1e4,1e5 = ['7.29e-02' '2.26e-02' '7.09e-03' '2.24e-03']
  [PASS] A1 atom weight at the lightest point -> 0 under refinement (a.c. cut, ~1/sqrt(nbins))
   free Gaussian: Var(phi)=0.9988 Var(phi^2)=1.9932  ratio=1.9979 (exact 2)
  [PASS] A2 Gaussian state obeys Wick (ratio=2): no bound-state atom
   => in EVERY Gaussian (Hartree) master field, at every beta incl. 2.5, g1/Sg = 0 exactly.
PART B  2D exact Singer W0 = exp(-lambda A(C)/2) : trF^2 = H conserved -> g1/Sg = 0
   SU(2) character Gram max|int chi_a chi_b - delta| = 2.22e-16
  [PASS] B1 SU(2) character orthogonality ; B2 trF^2 residues [0,0,0,0,~1e-16] -> g1/Sg=0
   generic 2D sources: holonomy chi_f -> g1/Sg=1.000 ; chi_f+chi_adj -> g1/Sg=0.500
  [PASS] B3 generic source: 2D W0 CAN give positive g1/Sg (trF^2 zero is a SELECTION RULE, not generic)
PART C  variational large-N ansatz : rigid-string breathing mode HAS binding
   kappa   rho*    m1=sqrt(E1-E0)   g1/Sg(O=V)   g1/Sg(O=rho)
   0.138   0.263   4.4956           0.9504       0.8949
   0.193   0.311   4.4137           0.9666       0.8986
  [PASS] C1 g1/Sg>0 STRICTLY (binding) ; C2 SOURCE-DEPENDENT (missing input) ;
         C3 rigid-string m1/sqrt(sigma)=2*pi*sqrt(2kappa) in band 3.3-3.9 (k=0.138-0.193)
PART D  the squeeze/ceiling bracket at beta=2.5 and what it does NOT give
   measure A (s=1:0.4,9:0.6) g1/Sg=0.400 ; measure B (s=1:0.001,2.5,6,12) g1/Sg=0.001 ; (M0,M1,M2) match
  [PASS] D1 equal moments ; D2 g1/Sg differs 400x => moments CANNOT bound g1/Sg (W2 NG)
# 10/10 checks PASS
```

Readings: **(A)** the Gaussian master field's bilinear source is a two-gluon bubble — the atom weight decays to
`0` under refinement and Wick (`ratio=2`) forbids a residue: `g₁/Σg=0` at every `β` `[T]`. **(B)** the 2D exact
`W₀` gives `0` for the conserved `tr F²` (orthogonality `2e−16`) but `1`/`1/2` for a generic source — a positive
control. **(C)** the rigid-string breathing truncation has binding with `g₁/Σg>0`, but source-dependent and
one-tube (upper estimate). **(D)** the moment-blindness no-go (`0.400` vs `0.001` at equal `(M₀,M₁,M₂)`).
**Witness, not proof; the theorems are §2–§5, the numerics illustrate the no-go, the exact anchors, and the
source ambiguity.**

---

## 8 · Ledger

| item | statement | label |
|---|---|---|
| Gaussian `W₀` construction | `W₀(C)=exp(½∮∮G)`, `G=1/(k²+μ²)`, `μ(β)` from tadpole — closes at every `β` | **[T]** |
| **Gaussian atom fraction** | bilinear `O=tr F²` ⇒ `⟨OO⟩_c=2G²` (bubble) ⇒ `dν` a.c. ⇒ `g₁/Σg=0` **∀β incl. 2.5** | **[T]** |
| Gaussian no-go | no self-consistent `μ(β)` (any `β`) makes a Gaussian `W₀` give `g₁/Σg>0` | **[NG]** |
| 2D exact `W₀` | `W₀=e^{−λA/2}`; `tr F²=H` conserved ⇒ `g₁/Σg=0`, all `N` | **[T]** |
| 2D generic source | `g₁/Σg=|o_{R₁}|²/Σ_{R≠0}|o_R|²>0`: holonomy `1`, `χ_f+χ_adj` `1/2` — a positive control | **[T]** |
| rigid-string `W₀` | binds: `m₁=2π√(2σκ)`, `m₁/√σ∈[3.3,3.9]` for `κ∈[0.138,0.193]` | **[C]** |
| rigid-string atom fraction | breathing truncation: `g₁/Σg=0.89–0.98` (`O=V`), `0.89–0.92` (`O=ρ`), `>0` ∀κ | **[V]** / **[C]** |
| rigid-string uncertainty | value **source-dependent** (`O=V` vs `O=ρ`) + one-tube (upper) ⇒ not a certified floor | **[C]**/**[O]** |
| source–tube coupling | `tr F²`↔tube functional **not fixed** by `W₀` — the missing input | **[O]** |
| `β=2.5` bracket | `0 ≤ g₁/Σg ≤ ceiling`; lower face from Gaussian/2D (exactly `0`), upper from squeeze + WG `χ_O` | **[T]** |
| WG ceiling | `χ_O=0.33–0.41` at `β=2.5` (gauge-invariant, non-self-averaging) bounds `g₁/Σg` **above** | **[V]** |
| moment blindness | `g₁/Σg` a point mass; equal `(M₀,M₁,M₂)` ⇒ `0.400` vs `0.001` — no moment functional bounds it | **[NG]** |
| **determination at `β=2.5`** | **no controlled truncation certifies a positive floor**; the floor `g₁/Σg≥w₀` **=** the non-Gaussian binding `= H1` | **[O]** |
| physical expectation | `g₁/Σg=O(0.2–0.5)` from `β=2.2` anchor (`f₁²_bare=7.14`) + trace anomaly + confinement | **[C]**/**[V-anchor]** |
| **Clay mass gap** | unbuilt | **[O] · 0%** |

---

## 9 · Sharpest line

**The master field the W2 note left unbuilt is built here — and it answers the floor question by *where it
splits*.** The **Gaussian/Hartree `W₀`** is an exact free massive functional at every `β`; its bilinear `tr F²`
source is a **two-gluon bubble**, so by Wick the KL measure is absolutely continuous and **`g₁/Σg=0` at every
`β`, `2.5` included** — the free-field acid test at finite coupling, and a **no-go for the whole Gaussian class**
`[T]/[NG]`. The **2D Singer `W₀`** returns `0` for the conserved `tr F²` but `1`/`1/2` for a generic source — a
**positive control** proving an exactly solvable `W₀` *can* carry a positive atom fraction, so the `tr F²` zero is
a selection rule, not a generic obstruction `[T]`. The **rigid-string `W₀`** has the binding the others lack —
`m₁=2π√(2σκ)` — and its breathing-mode truncation gives **`g₁/Σg>0`**, but the value is **source-dependent** (the
`tr F²`↔tube coupling is undefined) and the one-tube denominator omits the multi-tube continuum, so it is an
**upper estimate, not a certified floor** `[C]/[O]`. Therefore **at `β=2.5` no controlled truncation certifies a
positive floor**: the closed master fields give exactly `0`, the binding models are either outside the cluster
radius or source-ambiguous, the truncations **bracket** `0 ≤ g₁/Σg ≤ ceiling`, and the floor `g₁/Σg≥w₀>0` **is**
the non-Gaussian connected `1PI` `0⁺⁺` binding — `H1` in the `N`-coordinate — which no Gaussian `W₀` can hold.
The atom is intrinsically non-Gaussian; constructing a positive `g₁/Σg` from `W₀` and solving the binding are the
same problem. No QED. Clay 0%.

---

*Refs (internal).* `W2-MASTER-FIELD.md` (§1 ratio identity, §2 closed MF/BS equation, §4 λ-uniformity, §6
moment no-go), `W-WAVE-CHARTER.md` §1 W2/§2 acid battery; `U3-LARGE-N-BINDING.md` (Theorem U3-1, §3 2D
anchors); `W1-STRONG-ALIGN.md` §4 (`uK≫1` at `2.5`), `W3-SMALL-VOLUME.md` (exact small-`V` `w₁=0`);
`unconditional/CONSTRUCT-MASTERFIELD.md` (§1 closures (i) Gaussian→DSE wall, (iv) rigid string `m=2π√(2σκ)`;
§4 `κ=H1`), `PNT-LARGEN-LOOP.md` §1–§4 (`(MM∞)` closed; `W₀` 2D-only), `KAPPA-LARGEN.md`, `SLOT-DSE-FIXPOINT`,
`PROOF-EFFECTIVE-STRING.md`, `TRANSMUTE-STRING-SPECTRUM.md`, `ROOF-LARGEN-MASTER.md`; `R5-EXACT-ANCHORS.md`
§1.1 (2D `tr F²=H` atom `0`), §4 (binding is the discriminant); `T2-TUBE-OVERLAP.md` §8 (floor ⟺ isolation,
`w₀=O(0.2–0.5)`), `WG-INVARIANT-ESTIMATOR.md` §5–§6 (ceiling `χ_O=0.33–0.41`, squeeze), `Q2-LOOP-BOOTSTRAP.md`
§2,`docs/paper-v2/ch02-squeeze.tex` (`eq:squeeze`,`eq:hstar`). *External, standard, from the corpus (not
re-fetched):* 't Hooft NPB 72 (1974) 461; Makeenko–Migdal PLB 88 (1979) 135; Migdal Sov. Phys. JETP 42 (1976)
413; Singer (2D master field); Giles PRD 24 (1981) 2160; Källén (1952)/Lehmann (1954); Wick's theorem (Gaussian
measures); van Hove threshold. **Witness:** `research/v3/w7/w2b_master_w0.py` (`py -3.12 -I`, numpy,
deterministic, **10/10 PASS**, output `w2b_master_w0.out`). **No third-party downloads; no `[C]` upgraded to
`[T]`; no QED; Clay 0%.**
