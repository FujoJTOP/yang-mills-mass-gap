# R1-RESIDUE-TRANSPORT — attacking the atom, part R1: **what `Z_O` is, why `Z_O>0` is free and beside the point, and the exact missing step `Ω` (persistence = isolation `I₂` = `(CP)`), not `Z_O>0` and not a Ward identity**

> **2026-11-24 · Wave R1 — the multiplicative residue transport `f₁²(μ)=Z_O² f₁²(μ₀)` under the microscope.**
> Task: the atom target `f₁² = |⟨0|tr F²|G⟩|² = ν({m₁²}) ≥ r₀ > 0` with isolation `I₂`
> (`Λ_★ − m₁² ≥ δ > 0`). The corpus carries the P3 strong-coupling seed (`f₁²>0` at `β<1/48`,
> unconditional, `L`-uniform) to the physical head by the *multiplicative identity* `f₁²=Z_O² f₁²`.
> **Question:** what exactly is `Z_O`, can it be bounded below uniformly in `β` (especially at `β=2.5`),
> and can the strong-coupling residue be transported nonperturbatively? If conditionally, state the
> theorem + constants; else the exact no-go.
>
> **Verdict (no QED, Clay 0%). Three sharp answers, one of them a correction of emphasis.**
> **(1)** `Z_O` is the operator wave-function-renormalization of `O=tr F²`, fixed by the trace anomaly to
> `Z_O=(β_YM/2g)(μ₀)/(β_YM/2g)(μ)>0`; the identity is an **exact RG covariance within one continuum
> theory** `[T]`. **(2)** `Z_O` **is bounded below by a positive constant on every finite window**
> `[T]` — trivially, because `β_YM(g)<0` for all `g>0` in pure YM (no Banks–Zaks), so it is a ratio of
> two like-signed couplings; its **size** is scheme- and coupling-dependent and is a perturbative
> absurdity at `β_eff=1/48` (`g²=192`). **(3)** Bounding `Z_O` below **does not transport the atom**:
> the map is multiplicative and fixes `0`, so `Z_O²·0=0`. The missing datum is not `Z_O>0` and not a
> Ward identity (the Ward identity **is** `Z_O` and carries no lower bound on the amplitude); it is the
> **state-overlap / persistence `Ω(β₂,β₁)=|⟨G(β₂)|G(β₁)⟩_O|²`**, which equals `1` iff the same lightest
> `0⁺⁺` atom persists across the crossover — i.e. **uniform isolation `I₂`**, i.e. **`(CP)`**, i.e. the
> gap. `Ω` is an overlap, bounded only *from above* by Cauchy–Schwarz `[NG]` for any positivity/`ν≥0`
> class. Solvable cases (one-plaquette, 2D, strong-coupling) show the "transport constant" is **not
> universal** and sign-transport is **vacuous** at finite `(L,a)`. **No QED. Clay 0%.**
>
> **Reads.** `Q-WAVE-SYNTHESIS.md`, `Q6-H1-STATUS.md`; `PROOF-P3-UNCONDITIONAL.md`;
> `unconditional/INVENT-TRANSPORT-INEQ.md` (`Z_O` definition), `MS-RG-ODE.md` (the flow),
> `MS-FIXEDPOINT-RECURSION.md` (the signed defect `𝔪`), `OTHER-F1SQ-TRANSPORT.md`,
> `PROOF-CP-COMBINE.md`, `CERT-CP-TRANSPORT.md` (`[CP]=I₂=A3=H1`), `PROOF-COMPLEX-BETA-2D.md`
> (Kato: a carried path *is* the gap); `docs/paper-v2/ch02-squeeze.tex` (`eq:squeeze`, `eq:hstar`,
> `§schur`), `ch03-existence.tex` (`§dim`, `Z_O` transport); `lean/YMStatement/Reduction.lean`
> (`CP`, `CP_holds` 1 `sorry`, `transported_residue` 0 `sorry`).
> **Witness.** `research/v3/w7/r1_transport.py` (`py -3.12`, numpy/scipy, deterministic, one job).
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional · **[NG]** no-go for a
> named class · **[O]** exact missing step · **[V]** numerical witness here. No QED. Clay 0%.

---

## §0 · Verdict up front

| question | answer | label |
|---|---|---|
| What is `Z_O`? | the wave-function renormalization of the dim-4 scalar `O=tr F²`; the trace anomaly fixes `(β_YM/2g)O` RG-invariant ⇒ `Z_O=(β_YM/2g)(μ₀)/(β_YM/2g)(μ)>0`; `f₁²` runs `∝Z_O²` | **[T]** |
| Is the identity `f₁²(μ)=Z_O²f₁²(μ₀)` true? | yes, **as an RG covariance within one renormalized continuum theory, for the *same* state** | **[T]** |
| Can `Z_O` be bounded below uniformly in `β` (esp. `β=2.5`)? | **yes**, trivially: `Z_O²>0` is a square of a ratio of like-signed couplings; zero only at a fixed point (none in pure YM) or `μ→Λ` | **[T]** |
| Is the lower bound on `Z_O` the missing step? | **no**: it is free and beside the point; the map fixes `0` | **[T]** |
| Can the strong-coupling residue be transported nonperturbatively? | **only conditionally on `(CP)`** = uniform isolation `I₂` = the state-overlap `Ω=1` | **[C]** |
| What is the precise missing step? | `Ω(β₂,β₁)=|⟨G(β₂)|G(β₁)⟩_O|²=1` (persistence), equivalently uniform `I₂`; **not** `Z_O>0`, **not** a Ward identity | **[O]** |
| Exact no-go for the transport-alone programme | any bound on `Z_O` alone cannot bound `f₁²` below, because `f₁²=0 ⇒ Z_O²f₁²=0`; and `Ω` is pinned only from *above* by Cauchy–Schwarz | **[NG]** |
| Clay mass gap | unbuilt | **0%** |

**One line.** `Z_O>0` is trivial and was never the wall; the wall is `Ω` — that strong- and weak-coupling
sittings share the *same* lightest `0⁺⁺` atom — which is isolation `I₂`, i.e. `(CP)`, i.e. the gap. No QED.

---

## §1 · What `Z_O` is, exactly — the operator renormalization and the trace-anomaly covariance

### 1.1 The object

`O_S := ¼ tr F²` (any fixed normalization) is the dim-4, `0⁺⁺`, P-even scalar. Its connected (off-diagonal)
matrix element is multiplicatively renormalized:
```
   O_R(μ) = Z_O(μ)^{-1} O_B,        a_G(μ) := ⟨0|O_R(μ)|G⟩ = Z_O(μ)^{-1} a_G^{B},      (1.1)
```
The dim-≤4 partners of `O_S` are: the identity `1` (diagonal — it shifts the *vacuum expectation*
`⟨O⟩`, i.e. the gluon condensate, but contributes nothing to `⟨0|O|G⟩` for `G≠0`), the P-odd `F F̃`
(decouples from the `0⁺⁺` channel), and BRST-exact derivatives (vanish on physical states). So the
**off-diagonal** matrix element renormalizes multiplicatively with no mixing `[T]` (standard dim-4
operator mixing; Collins, *Renormalization*). This is what makes a *single* number `Z_O` well defined
for the residue `f₁²=|⟨0|O|G⟩|²`.

### 1.2 Why the anomaly fixes `Z_O` — and what "the same state" means

The trace anomaly is exact, `T^μ_μ = (β_YM/2g) O_S`, and `T^μ_μ` (the divergence of the dilatation
current) is RG-invariant. Hence `(β_YM/2g) O_S` is RG-invariant, so `Z_O ∝ β_YM/2g`:
```
   Z_O(μ,μ₀) = (β_YM/2g)(μ₀)/(β_YM/2g)(μ) > 0 ,      f₁²(μ) = Z_O(μ,μ₀)² f₁²(μ₀).     (1.2)
```
**The logical content of (1.2).** It is a statement about **one** renormalized continuum theory: the
matrix element of `O_R(μ)` between the **vacuum and a fixed physical state `|G⟩`**, evaluated at two RG
scales `μ₀,μ`. The state `|G⟩` does not depend on `μ` (it is a fixed vector in the physical Hilbert space;
the RG scale is a bookkeeping device for the operator). Therefore (1.2) is exact and carries **no**
dynamical information about *whether* `|G⟩` exists, is isolated, or is the lightest. That is all `[T]`.

**The correction of emphasis (§0).** This is precisely where the corpus's Lean `CP` bundles two very
different things (`lean/YMStatement/Reduction.lean` L164–166):
```
   def CP :=  0 < T.ZO   ∧   T.f1sq = T.ZO^2 * T.f1sqStrong
              └─ free ─┘    └──────── the hard equality ────────┘
```
The first conjunct `Z_O>0` is **[T]** (§2). The second conjunct is **the state identification** — that
`f₁²_strong` (the `β<1/48` object) **is** `Z_O²`-times the physical `f₁²` — and *that* is the `[O]`
(sorry) content. Reading `CP` as "`Z_O>0`" understates it; reading it as the equality is correct. This
note makes the split explicit.

---

## §2 · Is `Z_O` bounded below uniformly in `β`? — Yes, trivially; the size is not universal

### 2.1 The bound is free `[T]`

`Z_O(μ,μ₀) = (β_YM/2g)(μ₀)/(β_YM/2g)(μ)`. In pure Yang–Mills, `β_YM(g)<0` for every `g>0` (asymptotic
freedom, confirmed by the lattice β-function; **no Banks–Zaks zero**). Hence `β_YM/2g<0` at both ends,
the ratio is `>0`, and `Z_O²` is a square. The **only** zeros of `Z_O` are:
- a fixed point `β_YM(g_*)=0` — none exists in pure YM `[T-standard]`;
- `μ→Λ_QCD⁺` (deep IR, `g→∞`), where `f₁²→0` anyway.

So on any compact window `[μ_1,μ_2]⊂(Λ,∞)`, `Z_O ≥ c(μ_1,μ_2) > 0`. **`Z_O>0` and `Z_O` bounded below
are not difficulties.** Witness `[V]`: `r1_transport.py [C]` prints the tree-level value and the
`β_YM/2g` logic; `[A]` prints a concrete solvable-model analogue staying positive across the whole β-range.

### 2.2 The *size* is scheme/coupling-dependent — and meaningless at `β_eff` `[V]`

The numerical size of `Z_O` is **not** a universal datum:
- **tree level** `g²=2N/β` (N=2): `Z_O = g²(β_eff)/g²(2.5) = 192/1.6 = 120`, `Z_O² = 1.44×10^4`;
- **one loop** `1/g²=2b₀\ln(μ/Λ)`: `Z_O = \ln(\mu_{2.5}/Λ)/\ln(\mu_{eff}/Λ)`, a **log ratio**, a
  *different functional form*;
- **strong coupling** (§3): the residue ratio scales as a **power** `(β₂/β₁)^p`, again different.

All three are positive (so sign-transport is robust `[T]`), but **they disagree on the value**. The
tree-level `120` is computed from `g²(1/48)=192`, which is **non-perturbative** — the one-loop formula
that *defines* `Z_O` is invalid there. This is the honest reading: **`Z_O>0` is source-independent; the
magnitude is source-dependent and trustworthy only inside the scaling window.** The gap between them
is the crossover. Witness `[V]`: `r1_transport.py [B],[C]`.

---

## §3 · Solvable cases — the "transport constant" is not universal, and sign-transport is vacuous

### 3.1 One-plaquette / SU(2) character coefficients `[T]+[V]`

The one-plaquette model with action `S=β(1−cos(θ/2))`, `dU=π^{-1}\sin²(θ/2)dθ`, has the exact heat-kernel
(character) expansion `e^{β\cos(θ/2)}=Σ_j c_j(β)χ_j(θ)` with
```
   c_j(β) = 2(2j+1) I_{2j+1}(β)/β ,       (I_ν modified Bessel).        (3.1)
```
The natural "residue ratio" of the lightest excitation is `r(β):=c_{1/2}/c_0 = ⟨tr U⟩`. The witness
**verifies (3.1) against direct quadrature** (`Z₁=c₀` and `⟨tr U⟩=c_{1/2}/c_0` to 6 digits; `[A]`), then:

| `β` | 0.05 | 0.10 | 0.25 | 0.50 | 1.0 | 2.5 | 5.0 | 10 |
|---|---|---|---|---|---|---|---|---|
| `r=c_{1/2}/c_0` | .025 | .050 | .125 | .247 | .480 | 1.014 | 1.439 | 1.708 |
| `d\ln r/d\ln β` | — | .9995 | .9976 | .9889 | .9571 | .8157 | .5041 | .2479 |

Three exact facts `[T]+[V]`:
- `r(β)>0` for every `β>0` ⇒ **sign transport is vacuous** in every finite compact model (every residue
  is positive; the atom vs the cut is invisible to a sign).
- small `β`: `r∼β/2` ⇒ a **power law** `f₁²∼r²∼β²/4`; the local log-derivative `→1` at strong coupling,
  `→0` at large `β` — **not** a single `Z_O` power.
- `r(β₂)/r(β₁)` is a transcendental ratio, **not** the continuum anomaly `Z_O`.

### 3.2 Strong-coupling series `[T]+[V]`

At strong coupling `f₁²(β)=½c₁(β/8)^p(1+O(β))` (P3 leading tube `[T]`). The ratio between two couplings
is a **bare power** `(β₂/β₁)^p`, e.g. `(2.5/(1/48))²≈1.44×10⁴` at `p=2` (witness `[B]`). This is the
strong-coupling "transport constant"; it is **not** the anomaly's `β₂/β₁` (tree) or `[\ln]²` (1-loop).
The three constants disagree because their regimes do not overlap — **the crossover is exactly the
region none of them controls.**

### 3.3 2D (Migdal) `[T-standard]+[V]`

2D `SU(2)` YM is exactly solvable (Migdal). The transfer matrix in the representation basis is diagonal,
with the **heat-kernel weights `c_R(β)` of (3.1)** — the *same* Bessel object as §3.1. So **2D does not
supply a new transport constant; it reduces to the one-plaquette coefficients** `[T-standard]`. The
naive strong-coupling "mass" `−\ln(c_{1/2}/c_0)` is positive only for `β≲2` and turns negative at
`β=2.5` (witness `[F]`) — a concrete illustration that the strong-coupling expansion is *broken* at the
target coupling. (Caveat, as in the corpus: 2D is topological, no propagating gluon; it is a *solvable
laboratory* for the transport algebra, not a model of the 4D atom.)

**Reading of §3.** In every solvable case the "constant that carries `f₁²` between couplings" is
**scheme/coupling-dependent** (Bessel ratio / power / log-ratio), and sign-transport is **vacuous**
(`r>0` always, by compactness). Neither fact is the anomaly `Z_O`, and neither transports the atom as a
*floor*. The multiplicative law `Z_O²` is the *continuum scaling-window* statement, valid only where the
crossover is already assumed away.

---

## §4 · The exact decomposition: the transport carries the *value*, the missing factor is the *overlap*

### 4.1 The honest identity

Let `f₁²_{lat}(β)` be the lightest `0⁺⁺` atom of the *lattice* theory at coupling `β` (a residue of the
finite-`L` transfer matrix, which is always `>0` by compactness). To compare `β_eff<1/48` and `β_tgt=2.5`
we must identify the lattice theories as continuum regularizations and the states as the *same* physical
glueball. Inserting the two identifications:
```
   f₁²_{phys}(β_tgt) = Z_O²(μ_tgt,μ_eff) · Ω(β_tgt,β_eff) · f₁²_{lat}(β_eff) ,        (4.1)
   Ω(β_tgt,β_eff) := |⟨G(β_tgt)|G(β_eff)⟩_O|²   (O-weighted overlap of the two lightest atoms).  (4.2)
```
Here `Z_O²` is the exact §1 covariance (it acts on the operator), and `Ω` is the **state** factor — the
squared `O`-overlap between the lightest `0⁺⁺` eigenstate at the two couplings. The domain of (4.1)
requires both ends to be in the scaling window and the state to be the same. The literature identity
`f₁²=Z_O² f₁²` used by the corpus is the special case `Ω=1`.

### 4.2 What `Ω` is, and why it is exactly the missing step

`Ω` is `1` **iff** the strong-coupling "glueball" is the physical glueball, carried without mixing across
the crossover. If the lightest `0⁺⁺` at `β_eff` maps to a higher level, or dissolves into the continuum,
or crosses another level at `β_tgt`, then `Ω<1` and can be `0`. In KL coordinates:
- `Ω=1` ⇔ the lightest atom is the *same* atom, with the same quantum numbers, isolated on `[β_eff,β_tgt]`
  ⇒ **uniform isolation `I₂`** (`Λ_β−m₁²(β)≥δ>0`);
- `Ω=0` ⇔ the `O`-matrix element of the lightest state vanishes (an **overlap zero**, `a(β_tgt)=0`), or
  the atom has dissolved — `f₁²` touches `0`.

So **`Ω` is isolation `I₂`, is the single-`T=0`-phase statement `(CP)`, and is the gap.** The corpus's
own chain says exactly this (`OTHER-F1SQ-TRANSPORT`: transport *is* `(CP)`; `CERT-CP-TRANSPORT`: the
failing step is the `L→∞,a→0` survival = `I₂` = `H1`; `PROOF-COMPLEX-BETA-2D`: a Kato-carried path
carrying `f₁²` *is* the gap on `[1/48,2.5]`). This note sharpens the *name* of the missing factor to `Ω`,
separating it cleanly from `Z_O`.

---

## §5 · The no-go: no lower bound on `Z_O` can transport the atom `[NG]`

Three independent, elementary obstructions, each exact:

**(N1) `0` is an absorbing fixed point.** The map `f₁² ↦ Z_O² f₁²` is linear and vanishes at `0`, for
**every** `Z_O>0`. Hence `f₁²(β_eff)=0 ⇒ f₁²(β_tgt)=0`, and `f₁²(β_tgt)>0` carries **no** information
back-propagated. A "lower bound on `Z_O`" (`Z_O≥c>0`) multiplied by anything can be `0`. Witness `[D]`
(`r1_transport.py`): the two-level toy `H=[[a,v],[v,−a]]`, `O=σ_z` has residue
`f₁²=|\langle e_0|σ_z|e_1\rangle|²` that **vanishes at `v=0`** (selection rule) and is smooth/positive
off it; `Z_O²·f₁²→0` for any `Z_O`. The residue is *manufacturable-to-zero*; a multiplicative constant
cannot stop it.

**(N2) `Ω` is bounded from the wrong side.** By Cauchy–Schwarz / positivity, `|⟨G'|O|G⟩|² ≤
\langle G'|O²|G'\rangle\langle G|O²|G\rangle` — an **upper** bound on the overlap. Every `ν≥0`/RP
functional is Schur-concave (ch02 `§schur`), giving upper bounds on atom weights, never lower. So `Ω≥c>0`
is **not** in the positivity/RP/spectral class; it is a non-perturbative persistence statement, of the
same difficulty as the gap. `[NG]`.

**(N3) The Ward identity is `Z_O`, and carries no amplitude bound.** For the `0⁺⁺` glueball, the trace
anomaly read as a Ward identity gives
```
   ⟨0|T^μ_μ|G⟩ = m_G² f_G  =  (β_YM/2g) a_G ,    ⇒    a_G = (2g/β_YM) m_G² f_G ,        (5.1)
```
where `f_G` is the (dynamical, un-normalized) decay constant. Equation (5.1) **is** the source of the
multiplicative law (§1.2): `a_G ∝ (β_YM/2g)^{-1}`. It fixes the *scaling* and nothing else — it contains
`f_G`, a free matrix element. There is **no** Ward identity that lower-bounds `a_G` (equivalently `f₁²`),
because `f_G` and `m_G` are not bounded by the anomaly. So "transport via a Ward identity" is
**already used and already insufficient**: the Ward identity is `Z_O`, and `Z_O` is free.

**Conclusion of §5.** The transport-alone programme—"prove `Z_O>0`/bounded-below and ship the residue"—is
**refuted `[NG]`**: it is a covariance that fixes `0`, whose missing factor `Ω` is bounded on the wrong
side by every available positivity tool, and whose natural Ward-identity source carries no amplitude
bound. `Z_O>0` is necessary-shaped and trivially true; it is not sufficient.

---

## §6 · The conditional theorem and constants (state it honestly)

> **Theorem R1 (conditional residue transport) — `[T]` given `(CP)`.**
> Fix `N`, `T=0`, and a lattice/continuum scheme. Assume
> **(CP)**: (i) `β_eff<1/48` and the physical point `β_tgt` lie on one constant-physics line (a single
> confining phase), and (ii) the lightest `0⁺⁺` atom is Kato-isolated and carried throughout, i.e.
> `Ω(β_tgt,β_eff)=1` ⇔ uniform `Λ_β−m₁²(β)≥δ>0` on `[β_eff,β_tgt]`. Then
> ```
>   f₁²(β_tgt) = Z_O²(β_tgt,β_eff) · f₁²(β_eff) ≥ Z_O² · r₀ > 0 ,     Z_O=(β_YM/2g)(μ_eff)/(β_YM/2g)(μ_tgt)>0,
> ```
> with `r₀(β)` the P3 strong-coupling floor (`PROOF-P3-UNCONDITIONAL`, uniform in `L`) and
> `Z_O² ≥ c₀ > 0` an explicit positive constant (no Banks–Zaks zero). **Constants.** `r₀=½c₁(β/8)^p>0`
> (P3); `Z_O²` is the square of the anomaly-coefficient ratio (=`[\ln]_{ratio}²` at one loop, `>1`
> toward the UV), `(L,a)`- and `V`-uniform (a pure ratio of running couplings, never `e^{−cV}`).
> **This is a `[T]`-given-`(CP)` statement; `(CP)` is NOT proved and is exactly the gap.**

> **Corollary (what is *not* available).** `[NG]` No unconditional `f₁²(β_tgt)>0` follows from P3 plus the
> multiplicative identity alone, because `Ω` can be `0` (N1) and is bounded only from above (N2). The
> unconditional residue at `β=2.5` is **not** obtained.

**Constants, explicitly and honestly.** `Z_O²` is well-defined and positive; its numerical magnitude is
regime-dependent (§2.2, §3) and trustworthy only in the scaling window. `r₀` is the P3 floor. Their
product `Z_O² r₀` is a **conditional** lower bound, meaningful only once `(CP)` (hence the scaling
window) holds. Quoting a number `Z_O² r₀` at `β=2.5` would require `(CP)` — i.e. the gap.

---

## §7 · Acid tests `[T]+[V]`

1. **Free massless field.** `f₁²=0` (no atom; `:φ²:` has `a≡0` by `Z₂`). `Z_O²·0=0`. The map
   manufactures **no** atom — and *that* is why it is acid-clean **and** why it cannot certify the gap.
   `[T]`, witness `[E]`.
2. **`:φ²:` massive free composite.** Gapped (`m₁²=4m²>0`) but `f₁²=0` (a cut, not an atom):
   `Z_O²·0=0`. Consistent. `[T]`.
3. **Multi-atom / cut.** `Z_O²` multiplies every spectral weight alike, so the *shape* (atom vs cut) is
   preserved, not created: `Ω` is a property of the support, invisible to `Z_O`. `[T]`.
4. **One-plaquette / 2D / strong-coupling finite models.** `r(β)>0` for all `β` (§3) — sign-transport is
   vacuous, the discriminator (atom vs cut) is not addressed. `[V]`.

---

## §8 · The precise missing step, and what would close it

**Missing step (exact).** A non-perturbative lower bound
```
   Ω(β_tgt,β_eff) = |⟨G(β_tgt)|O|G(β_eff)⟩|² ≥ ω₀ > 0 ,        (P-Ω)
```
equivalently **uniform isolation** `I₂`: `Λ_β − m₁²(β) ≥ δ>0` along `[β_eff,β_tgt]`, equivalently the
lightest `0⁺⁺` atom persists without dissolution or level-crossing. This is **not** `Z_O>0`, **not** a
Ward identity (the Ward identity is `Z_O`, §5-N3), and **not** a positivity/moment inequality (which
gives only upper bounds on overlaps, §5-N2). It is the **gap** in residue costume: uniform `I₂`
(`Λ_★−m₁²≥δ`) is, with `Λ_★≈(1.5m₁)²` the excited-`0⁺⁺` threshold, equivalent to `m₁≥μ>0`
(`PROOF-CP-COMBINE §3`, `OTHER-F1SQ-TRANSPORT`).

**What would close it (structural).** Any of: (a) a continuum RG map that is a *contraction in the
state space* carrying the lightest `0⁺⁺` eigenvector with it (a Kato-stable "same-atom" flow with a
uniform spectral-isolation margin); (b) a direct lattice measurement of `Ω` with certified isolation
`I₂` on the `0⁺⁺` channel at `β=2.5` (a decisive-numerics target, not available); (c) a proof that the
strong-coupling and physical `0⁺⁺` atoms are connected by an adiabatic no-crossing path — which is Kato
analyticity on an interval that *is* `I₂` (`PROOF-COMPLEX-BETA-2D`). All three are the gap. **The
transport (`Z_O`) is free; `(P-Ω)` is the work, and it is unmoved.**

---

## §9 · Status ledger

| item | statement | label |
|---|---|---|
| `Z_O` = wave-function renorm. of `O=tr F²`, dim-4 no-mixing (off-diagonal) | §1.1; Collins | **[T]** |
| `(β_YM/2g)O` RG-invariant ⇒ `Z_O=(β_YM/2g)(μ₀)/(β_YM/2g)(μ)>0`; `f₁²∝Z_O²` | (1.2); trace anomaly | **[T]** |
| identity valid only **within one continuum theory, same state** | §1.2 | **[T]** |
| **`Z_O` bounded below by `c₀>0`** on any window (no Banks–Zaks zero) | §2.1; `[C]§` logic | **[T]** |
| `Z_O` *size* scheme/coupling-dependent; tree=120, 1-loop=log, SC=power | §2.2, §3; `[B],[C]` | **[T]+[V]** |
| one-plaquette `c_j=2(2j+1)I_{2j+1}/β`, `r=c_{1/2}/c_0>0 ∀β` | (3.1); `[A]` verifies vs quadrature | **[T]+[V]** |
| sign transport **vacuous** at finite `(L,a)` (every residue `>0`) | §3.1; CERT-CP (1) | **[T]+[V]** |
| strong-coupling transport is a **power** `(β₂/β₁)^p`, ≠ anomaly | §3.2; `[B]` | **[T]+[V]** |
| 2D (Migdal) weights = same `c_j(β)`; no new constant; SC mass breaks at `β=2.5` | §3.3; `[F]` | **[T-standard]+[V]** |
| decomposition `f₁²(tgt)=Z_O²·Ω·f₁²(eff)`, `Ω=|<G_tgt|O|G_eff>|²` | (4.1); `Ω=1`⇔same atom | **[T]** |
| `Ω = isolation I₂ = (CP) = gap` | §4.2; OTHER-F1SQ, CERT-CP, COMPLEX-BETA-2D | **[T]** (identity of objects) |
| **no lower bound on `Z_O` transports the atom**: `0` absorbing, `Ω` upper-bounded | §5 N1,N2; `[D]` | **[NG]** |
| Ward identity = `Z_O`, no amplitude bound (`f_G` free) | §5 (5.1) N3 | **[T]/[NG]** |
| conditional theorem `f₁²(tgt)≥Z_O²r₀>0` **given `(CP)`** | §6; Lean `transported_residue` 0 sorry | **[C] given (CP)** |
| unconditional `f₁²(2.5)>0` | **not obtained** | **[O]** |
| **precise missing step** `(P-Ω)`: `Ω≥ω₀>0` ⇔ uniform `I₂` ⇔ `m₁≥μ` | §8 | **[O]** |
| Clay mass gap | unbuilt | **0%** |

---

## §10 · Sharpest line

**`Z_O` is the wave-function renormalization of `tr F²` fixed by the trace anomaly, the identity
`f₁²(μ)=Z_O²f₁²(μ₀)` is an exact RG covariance within one continuum theory, and `Z_O>0` — indeed
`Z_O` bounded below by a positive constant — is trivially true in pure YM (no Banks–Zaks zero), so it was
never the wall.** The wall is the factor the covariance does **not** carry: the state-overlap
`Ω(β_tgt,β_eff)=|⟨G(β_tgt)|O|G(β_eff)⟩|²`, which equals `1` iff the strong-coupling lightest `0⁺⁺` atom
**is** the physical one — i.e. **uniform isolation `I₂`**, i.e. `(CP)`, i.e. the gap. The multiplicative
map is a homomorphism on the *weight* that fixes `0`, so no bound on `Z_O` (N1) can lift a dissolved atom;
and `Ω` is bounded only **from above** by positivity/Cauchy–Schwarz (N2), so no `ν≥0`/RP/Schur tool can
supply it; and the Ward identity that *is* the source of `Z_O` carries the free decay constant `f_G`
(N3) and no amplitude bound. Solvable cases (one-plaquette, 2D Migdal, strong-coupling series) show the
"transport constant" is **scheme/coupling-dependent** and sign-transport is **vacuous** at finite
`(L,a)`. **So: transport works conditionally — `f₁²(tgt)≥Z_O²r₀>0` given `(CP)` — and the exact missing
step is `(P-Ω): Ω≥ω₀>0`, equivalently uniform `I₂`, equivalently the gap; it is not `Z_O>0`, not
isolation-in-Z_O's-clothing, and not a Ward identity. No QED. Clay 0%.** The irreducible target is
untouched: the atom persists or it doesn't, and no multiplicative constant decides that.

---

*Refs.* `INVENT-TRANSPORT-INEQ.md` (`Z_O`, `γ_O`, `0→0` acid); `MS-RG-ODE.md` (`d\ln f₁²/d\lnμ=2γ_O>0`,
`m₁` RG-invariant, ratios RG-invariant); `MS-FIXEDPOINT-RECURSION.md` (signed defect `𝔪`, `c_eff=c+𝔪`,
floor⇔`I₂`); `OTHER-F1SQ-TRANSPORT.md` (transport is `(CP)`; failure is `I₂` not single-phase);
`PROOF-CP-COMBINE.md` (`I₂`: `Λ−m₁²≥κ`⟺`m₁≥μ`; `2m₁` vs `1.5m₁` threshold; overlap-zero `a(2.5)=0`);
`CERT-CP-TRANSPORT.md` (finite-`(L,a)` vacuous; flow downstream; renorm local ⇒ blind to `I₂`);
`PROOF-COMPLEX-BETA-2D.md` (Kato: a carried path *is* the gap on `[1/48,2.5]`); `SLOT-RICCATI-FLOW.md`
(`0` absorbing; linear flow sign-preserving, cannot manufacture); `PROOF-P3-UNCONDITIONAL.md` (`r₀>0`,
`β<1/48`, `L`-uniform; explicitly retracts `I₂`); `docs/paper-v2/ch02-squeeze.tex` (`eq:squeeze`,
`eq:hstar`, `§schur` Schur-concavity ⇒ upper-only); `ch03-existence.tex` (`§dim`, `Z_O²` dimension); 
`lean/YMStatement/Reduction.lean` (`CP` L164–166 splits `Z_O>0` from the equality; `CP_holds` 1 `sorry`;
`transported_residue` 0 `sorry`, conditional). **External (standard, cited from corpus, not re-fetched):**
Collins, *Renormalization* (dim-4 operator mixing, anomalous dimension of `tr F²`); Nielsen, NPB **120**
(1977) 212; Adler–Bardeen–Jackiw trace anomaly; Källén (1952)/Lehmann (1954) KL representation; Kato,
*Perturbation Theory for Linear Operators* (1976) Ch. VII §3 (isolated simple eigenvalue, analytic
transport); Osterwalder–Seiler, Ann. Phys. **110** (1978) 440 (RP transfer matrix); Migdal (1975) (2D
solvability); Shen–Zhu–Zhu, CMP **400** (2023) 805 (strong-coupling gap, `1/48`). **No third-party sources
downloaded; one witness run (`w7/r1_transport.py`, `py -3.12`, deterministic). No `[C]` upgraded to `[T]`;
no QED; Clay 0%.**
