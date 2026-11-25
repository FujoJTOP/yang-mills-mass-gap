# T5-CENTER-VORTEX — the centre-vortex bound on **binding**: the vortex picture is rigorous in the **charged** sector (Tomboulis–Yaffe `σ>0`) and **exactly blind** to the neutral `tr F²` atom; on the isolate ratio it gives a heuristic, not an inequality; the exact missing estimate is a *sign-cancelling positive regrouping of the dual-surface weight in the N-ality-0 sector* at β≈2.5 (≡ `(E-I2)`)

> **2026-11-24 · Wave T5 — centre-vortex species against the binding criterion.**
> The atom `f₁²=|⟨0|tr F²|G⟩|²` is a **binding** phenomenon; the minimal connected `tr F²`
> polymer is one four-walled flux tube (`P3-FOURWALL`, `[T]`). Centre vortices are the corpus's
> **confinement** picture (`PROOF-CENTER-VORTEX`, `probe-center/REPORT`). Question: can the tube
> overlap or the isolation ratio `ρ=m_{0⁺⁺*}/m_{0⁺⁺}` be bounded through the vortex representation?
> **Reads.** `R-WAVE-SYNTHESIS.md`, `R5-EXACT-ANCHORS.md`, `R2-ISOLATION.md`,
> `docs/paper-v2/ch04-gap.tex`; `unconditional/PROOF-GAUSS-INVARIANT-PF.md`,
> `PROOF-CENTER-VORTEX.md`, `PROOF-RANDOM-CLUSTER.md`, `PROOF-CONFINEMENT-RANDOM-SURFACE.md`,
> `P3-FOURWALL-DERIVATION.md`; `v3/P3-POLYAKOV-CENTER.md`; `research/probe-center/REPORT.md`.
> **Witness.** `research/v3/w7/t5_center_vortex.py` (`py -3.12 -I`, numpy only, deterministic, ~1 s;
> output quoted §8). No third-party downloads.
> **Labels (charter §2).** **[T]** true theorem / exact identity · **[C]** conditional ·
> **[NG]** no-go for a named class · **[O]** exactly the missing step · **[V]** numerical witness here.
> **No QED. Clay 0%.**

---

## 0 · Verdict (up front; no QED, Clay 0%)

**The centre-vortex picture does *not* give a rigorous bound on binding, and the failure is a
theorem, not an accident.** Three statements, in order of strength:

1. **[T] Exact neutral blindness.** `tr F²` is a colour singlet, **N-ality 0**. A centre vortex
   multiplies a Wilson loop in rep `R` by `χ_R(z)/d_R`, `z∈Z(N)`, and for **every** N-ality-0
   representation this factor is **exactly 1** (`χ_f`-ratio `= ω` vs adjoint ratio `= 1`; witness
   `[A]`, `max|neutral−1| = 4.4×10⁻¹⁶`). The Z(N) *phase* — the entire content of "centre dominance"
   and of vortex condensation — carries **zero** information into the neutral local sector where the
   atom lives. This is the group-theoretic core of `probe-center/REPORT` and `P3-POLYAKOV-CENTER`.

2. **[T]+[V] `σ>0 ⇏ atom`, theorem-by-example.** Confinement and the `tr F²` atom are *different
   objects*, and 2D pure Yang–Mills proves it in one exactly-solvable model: it **confines**
   (`σ_fund = g²C₂/2 > 0`, exact Migdal) **and** has `f₁²(tr F²) = 0` exactly (character
   orthogonality `⟨χ₀|χ_R⟩=δ_{0R}`, `tr F² = tr E² = H` conserved; witness `[B]`, `R5 §1.1`). So the
   only rigorous vortex route (`vortex condensation ⇒ σ>0`, Tomboulis–Yaffe) lands in the **1-form
   channel**, and `σ>0` is *provably* compatible with `f₁²=0`. The vortex route to the atom is dead.

3. **[NG] No neutral-sector percolation measure.** The only route that *could* bound the tube overlap
   from the vortex picture is a **positive** dual-surface (vortex) measure whose connection
   probability lower-bounds the tube amplitude. For `Z₂/Z_N/U(1)` this exists (Brydges–Giffen–Durhuus–
   Fröhlich positive weights; FK/FKG; `PROOF-RANDOM-CLUSTER`, `[T]`) — but it is **abelian only**. For
   `SU(N≥2)` in `d=4` the surface weights are the **Weingarten function, which is signed** (for SU(2)
   the sign is the Euler characteristic; for SU(N) nearly half are negative). There is **no positive
   Gibbs/FKG measure** at β≈2.5, hence no percolation inequality to state — a genuine no-go, not a gap.

**On the isolation ratio `ρ`.** Even granting a positive tube amplitude, `ρ` needs the *second*
closed-surface mass — an **order-2** datum. The vortex picture supplies no inequality for the ratio of
the two lightest closed-surface masses; that is the same order-2 (`TP₂`/second-cone) wall R2 isolated,
in surface language (`R2 §3.3`). The vortex representation moves the wall, it does not lower it.

**One line.** *Centre vortices are rigorous exactly one way — Tomboulis–Yaffe turn "vortices condense"
into `σ>0` — and that is the charged channel, provably blind to the neutral `tr F²` atom (2D: `σ>0`,
`f₁²=0`). The neutral sector has no positive vortex measure in non-abelian 4D (Weingarten signs), so the
tube overlap has no vortex lower bound and `ρ` has no vortex inequality. The exact missing estimate is a
sign-cancelling positive regrouping of the dual-surface weight in the N-ality-0 sector at β≈2.5 — the
same object as `(E-I2)`. No QED. Clay 0%.*

---

## 1 · The centre-vortex expansion of the `tr F²` correlator

We must keep apart **two things the phrase "vortex picture" conflates** (`P3-POLYAKOV-CENTER §1`,
`PROBE-CENTER §1–4`):

- **(a) the exact dual-surface expansion** — the character / random-surface representation of the full
  measure. This preserves *all* information; it is the polymer expansion in surface language.
- **(b) the centre-projected vortex condensate** — maximal-centre-gauge `U=zV`, `z∈Z(N)`, keeping only
  the centre part `z`. This is the *reduced mechanism* ("centre dominance", vortex percolation).

The two have completely different status, and conflating them is the source of the confusion the
corpus warned about ("vortices exist" ≠ "vortices carry the neutral sign", `probe-center/REPORT §0`).

### 1.1 The three levels of the expansion `[T]`

For `SU(N)` Wilson lattice YM, `S=βΣ_p(1−N⁻¹Re tr U_p)`, expand
`e^{(β/N)Re tr U_p}=Σ_r c_r(β)χ_r(U_p)`. Summing link-by-link, the connected `0⁺⁺` correlator becomes

```
  C(t) := ⟨tr F²(x,0) tr F²(y,t)⟩_c = Σ_{S : ∂S = P_x ∪ P_0}  W_v(S) ,          (1.1)
```

where `P_x, P_0` are the two source/sink plaquette loops, `S` runs over embedded **dual surfaces**
(the vortex worldsheets / character surfaces) whose boundary is those two loops, and `W_v(S)` is the
product of `c_r` and **Weingarten** recoupling weights of the surface tiling. This is the
Cao–Park–Sheffield exact representation (`arXiv:2307.06790`) restricted to two local insertions;
`W_v(S)` is **signed** for `SU(N)` (`PROOF-RANDOM-CLUSTER`). *(1.1) is [T] as a representation.*

- **Level (a).** (1.1) is exact and equals the polymer expansion. The Källén–Lehmann residue is
  `f₁² = |⟨Ω|tr F²|G⟩|² = |amplitude of the minimal spanning surface|²`; the lightest `G` is the
  lightest closed surface (the "vortex/flux glueball"). *This is the honest vortex expansion of the
  `tr F²` correlator.*

### 1.2 Which vortex configurations contribute to the minimal tube `[T]`

By `P3-FOURWALL §1`, the minimal connected `tr F²` polymer is the **four-walled straight tube**:
a cylinder of cross-section the plaquette boundary, four temporal plaquettes per slice, `|P|=4t/a`;
its (unique) leading class carries amplitude `a₀(β)= (2N)⁻²·[positive endpoint recouplings]·[finite
transverse sum] > 0`, and it is the world-line of the quenched plaquette flux. In the dual picture
(1.1) the configurations that contribute to *this* tube are exactly the **minimal spanning surfaces**
— surfaces of area `4·N_t` in the (dual) plaquette grading that coincide with the straight cylinder on
each time slice. Geometrically the minimal tube and the minimal vortex surface *coincide* (both are
cylinders over the plaquette); they are the **electric** and **magnetic** duals of one another. For the
**neutral** `tr F²` tube, both carry **N-ality 0**: the centre flux that would weight a vortex is
`z⁰ = 1`. So the centre charge contributes **nothing** to the minimal tube (§2).

### 1.3 Level (b): the centre projection is neutral-blind `[T]`

In (b), `U=zV`. The `tr F²` correlator factorizes into a centre-phase factor times a `V`-part
correlator; by §2 the centre-phase factor is **1**, so the entire neutral content sits in the `V`
part — precisely what centre projection **discards**. This is why "centre dominance" (a statement
about the *fundamental* Wilson loop, `probe-center/REPORT §4`: `χ_proj/χ_full≈1.5–1.7`) reproduces
the charged string tension and is **identically trivial in the adjoint** (`max|W_A−1|=0`).

---

## 2 · The exact obstruction: neutral blindness `[T]`

> **Lemma CV-1 (centre blindness of the N-ality-0 local sector) `[T]`.**
> For `z∈Z(N)` and any rep `R` of N-ality `q`, `χ_R(z·I)/d_R = ω^{q}`, `ω=e^{2πi/N}`. Hence for
> `q=0` — the adjoint, and **every colour singlet local operator** including `tr F²` — the phase is
> **exactly 1**. All centre-vortex phases act as the identity on `tr F²`.

*Proof.* A centre element `zI` is central, so it acts as the scalar `z^q = ω^q` on every weight of `R`;
thus `χ_R(zI)=d_R ω^q`, `d_R` cancels. ∎

*Scope note (why `tr F²` is N-ality 0).* `tr F²` is built from adjoint-valued `F^a_{μν}` and contracted
to a colour singlet; the **adjoint has N-ality 0**, so `χ_adj(zI)/d_adj = 1`. Equivalently, and more
invariantly, a **local** gauge-invariant operator carries **no 1-form centre charge** — a point has no
lattice linking with a vortex worldsheet (`PROOF-CENTER-VORTEX §4`: "the symmetry is provably blind to
the local `tr F²`"). (A *bare* lattice plaquette `tr U_p` is a small fundamental loop and does carry a
centre phase order-by-order in `a`; but the residue `f₁²` and the connected correlator of the
renormalised local operator are the neutral, physical objects, and the `a→0` linking probability of the
shrinking loop with a vortex vanishes.) The group computation `[A]` is exactly the adjoint statement.

Witness `[A]` (exact):
```
  SU(2) -I:  fundamental χ/d = -1.000000 ;  adjoint χ/d = +1.000000 ;  j=2 χ/d = +1.000000
  SU(3) ω^k: fund phase = ω^k (k=1: -0.5+0.866i) ;  adjoint(N-ality 0) phase = +1.000000 (all k)
  Z(N):      q=0 phase = +1 exactly ;  q=±1 phase = ω^{±1}
  max |neutral phase − 1| over SU(2),SU(3) adjoints = 4.4e-16
```

**Consequence.** The one channel by which a vortex could carry a sign into a real observable — the
`Z(N)` phase — is *exactly the identity* on the neutral sector. So:

- The signed/absolute vortex activity ratio `z_signed/z_abs` is the constant `−1/2` in the **charged**
  (fundamental) sector and **`+1`** in the **neutral** sector, for every cluster, every size
  (`probe-center/REPORT §3`, D1/D2). There is no size-dependent cancellation to build a bound on.
- "Vortex condensation", "centre dominance", "percolation" are all **phase/charged-sector** statements.
  They cannot manufacture a neutral-sector sign, exactly as `probe-center` concluded: *a mechanism that
  is exactly neutral-blind cannot supply the neutral sign cancellation the atom needs.*

This is why the vortex picture is *not merely unproved* for the atom — its central object (the centre
phase) is *provably trivial* on `tr F²`.

---

## 3 · The rigorous inequality that IS available — and why it cannot reach the atom `[T]`

The one rigorous centre-vortex inequality in the corpus is **Tomboulis–Yaffe** (`CMP 100 (1985) 313`,
`PROOF-CENTER-VORTEX §2`):

```
  ⟨W_fund(C)⟩ ≤ 2 { ½(1 − Z⁻_Λ/Z_Λ) }^{ A_C /(L₁L₂) } ,     Z⁻/Z = vortex-twisted partition fn.   (3.1)
```

Granting a vortex free-energy density of the condensing sign, (3.1) forces `⟨W_fund⟩` area-law, i.e.
`σ ≥ c·ρ̂ > 0`. This is **`[T]`, but**:

- **(i) Charged sector only.** `W_fund` is N-ality 1; (3.1) says nothing about `tr F²`. By Lemma CV-1
  the vortex free energy `Z⁻/Z` is a centre-sector (charged) object: it enters the neutral correlator
  as `1`.
- **(ii) `σ>0 ⇏ atom, exactly.** 2D pure YM is the decisive counterexample: it **confines**
  (`σ_fund=g²C₂/2>0`, exact Migdal, every β) and has `f₁²(tr F²)=0` exactly. Witness `[B]`:
  ```
     Haar orthogonality ∫χ_j dU = 1, 0, 3e-17, 3e-17, -7e-18  (δ_{j0})
     => ⟨0|tr F²|R⟩ = C₂(R)·δ_{0R} = 0  for all R, all t   => no m>0 atom
     β=1.0: σ_fund=1.50000>0 ;  β=2.5: 0.60000>0 ;  β=10: 0.15000>0   (CONFINES)
  ```
  So even a *proof* of `σ(2.5)>0` (still `[O]`, `PROOF-CENTER-VORTEX §3`) yields only the 1-form/flux
  mass `m₁=√σ`, **not** the `0⁺⁺` atom. The bridge `m_{0⁺⁺}=c√σ` (`c≈3.5–4`, Morningstar–Peardon) is a
  *measured* transmutation identity, not a theorem (`PROOF-CONFINEMENT-RANDOM-SURFACE §3`).

**Reading.** The centre-vortex picture is rigorous exactly one way (`criterion ⇒ σ>0`), and that way is
provably orthogonal to the atom (§2, §3(ii)). The "vortex `⇒` confinement `⇒` (measured) `m_{0⁺⁺}`"
chain has a **measured** link (`m∝√σ`) and a **trivial-on-neutral** link (Lemma CV-1): neither is a
binding bound.

---

## 4 · Attempted lower bound on the overlap and on `ρ` — where it collapses `[O]`/`[NG]`

We *attempt* the bound honestly, in the dual representation (1.1), and exhibit the exact failure.

### 4.1 The overlap `f₁²`

From (1.1), the residue is the modulus-square of the minimal-surface amplitude in the dual expansion:

```
  f₁² = | ⟨Ω|tr F²|G⟩ |² ,   ⟨Ω|tr F²|G⟩ = Σ_{S : ∂S = P_x∪P_0} W_v(S) .                 (4.1)
```

**Attempted lower bound.** Isolate the minimal tube `S★` (straight four-wall cylinder, §1.2):
```
  f₁² = | W_v(S★) + Σ_{S≠S★} W_v(S) |²  ≥  ( |W_v(S★)| − Σ_{S≠S★}|W_v(S)| )²_+ .            (4.2)
```
(4.2) is a genuine lower bound **iff** the absolute remainder is smaller than the leading amplitude,
`Σ_{S≠S★}|W_v(S)| < |W_v(S★)|`.

- **Strong coupling (`β<1/48`): it holds, and is exactly P3 `[T]`.** The cost grading
  (`P3-FOURWALL §3`) gives `Σ_{S≠S★}|W_v(S)| ≤ O(u)|W_v(S★)|` for bounded-excess surfaces, with
  `u(β)=c_f/N→0`. But the *decoration* sum is `t`-linear before resummation, so making (4.2) `t`-uniform
  is precisely the **dressed-mass/amplitude resummation `[R]`** of `P3-FOURWALL §3.2`: existence of
  the resummed weight, its `L`-uniformity, and a *lower* bound `a₀(1−O(u))>0` from a **signed**
  remainder (`|1+r|≤1+O(u)` does **not** give `1+r≥1−O(u)` without controlling the sign of the
  remainder). So even the positive case is `[C]` conditional on `[R]` — and it lives at `β<1/48`, far
  from 2.5.
- **β≈2.5, non-abelian: (4.2) collapses `[NG]`.** `W_v(S)` is the **Weingarten** weight — signed
  (`PROOF-RANDOM-CLUSTER`). For `SU(N)` nearly half the surface weights are negative, so
  `|Σ W_v|` can be arbitrarily smaller than `Σ|W_v|`; the bracket `|W_v(S★)|−Σ|W_v|` need not be
  positive, and (4.2) gives only `f₁² ≥ 0`. There is **no** positive Gibbs/FKG measure to repair this:
  the sign is topological (Euler characteristic for SU(2)), not a small-activity accident.

**The positivity that saves strong coupling is NOT vortex percolation.** In P3 the amplitude
`a₀>0` comes from the SU(N) **recoupling** `f⊗f̄→1` coefficient `1/N>0` (`P3-FOURWALL §2.1`), an
*algebraic* positivity — witness `[C]`: `u(β)=I₂/I₁>0` and `1/N=1/2>0`. This positivity is *independent
of the vortex phase* and would hold even if vortices were removed. It is not a statement about centre
flux; it is the recoupling algebra of the tube endpoints.

### 4.2 The isolation ratio `ρ`

Even granting `f₁²>0`, `ρ=m_{0⁺⁺*}/m_{0⁺⁺}` is the ratio of the two lightest **closed-surface**
(vortex-glueball) masses. The vortex representation writes them as the low-lying spectrum of the
closed-surface transfer operator — but the *ratio* is a **second-eigenvalue** datum. No vortex-based
inequality forces `ρ>1`: doing so needs an **order-2 (total-positivity / second-cone)** structure on the
surface ensemble (`R2 §3.3`, `PROOF-GAUSS-INVARIANT-PF`), which the orderless gauge quotient lacks.
The vortex picture therefore moves the order-2 wall from `𝒜/𝒢` to "closed surfaces" without lowering it.

- **Abelian (`Z(N)` 4D, 3D):** `ρ` is *computable* (positive dual weights; §5) but still requires
  solving the dual statistical model — no closed-form, no free inequality.
- **Non-abelian 4D β≈2.5:** `ρ` is `[O]` — the same `(E-I2)` of R2.

---

## 5 · Test battery: `Z(N)` / 3D / 2D / free field / strong coupling

Each test pins where the vortex picture is rigorous and where it is neutral-empty.

### 5.1 4D `Z(N)` lattice gauge theory — the vortex picture is **exact**, abelian `[T-standard]`+`[V]`

`Z₂` gauge in `d=4` is **Wegner-dual** to the 4D Ising model; `Z_N` to the `Z_N` clock/spin model. In
`Z(N)` the only content is the centre itself, so **the vortex picture is not a truncation — it is the
exact dual model** (no `V` part). Its surface weights are **positive** (FK/Edwards–Sokal), so the
gap ⟺ percolation criterion is a **theorem** (`PROOF-RANDOM-CLUSTER`). Hence in the confined phase the
`Z(N)` "`tr F²`" (plaquette) correlator is rigorously gapped and `f₁²>0` — *but*: (i) this is the
**abelian** luxury of positive weights; `Z₂`/`Z_N` have no N-ality-0 non-trivial rep, so the plaquette
*is* charged and blindness never bites; (ii) the *value* and `ρ` are not closed-form (still a 4D
statistical model). **Lesson:** the vortex picture is exact **iff** the dual weights are positive —
which is exactly the property `SU(N)` lacks (§5.5).

### 5.2 3D `Z(N)` / compact `U(1)` `[T-standard]`

3D `Z(N)` gauge theory is dual to the `Z(N)` clock model; compact `U(1)₃` confines by **monopole
condensation** (Polyakov), and the monopole world-line **is** the centre-vortex line. Here the vortex
picture = monopole picture and the mass gap **is** monopole/vortex-controlled (semiclassical,
computable). The `tr F²` (`F²` in 3D) correlator is a neutral composite but in this abelian 3D setting
the dual expansion is positive and the gap follows. **Lesson:** the cleanest case where "vortex `⇒`
mass gap" is *true and computable* — and it is abelian and 3D.

### 5.3 2D pure YM — the **acid/orthogonality** counterexample `[T]`+`[V]`

Exactly solvable (Migdal). It **confines** (`σ_fund=g²C₂/2>0` for all β) **and** has `f₁²(tr F²)=0`
exactly (`tr F²=tr E²=H` conserved; character orthogonality). Witness `[B]`. This is the single
strongest test: a theory in which the vortex/confinement mechanism is *fully operative* (area law) yet
the `tr F²` atom is *exactly zero*. **Any vortex-based binding bound that returns `f₁²>0` or `ρ>1` in
2D is refuted.** It is the `tr F²`-analogue of the free-field acid test (`R5 §1.1`).

### 5.4 Free massless field — acid test `[T]`+`[V]`

No vortices (trivial/flat measure); `tr F²` is a two-particle cut, `f₁²=0` (witness `[D]`). A vortex
bound must return `0`. It does, trivially (the dual surface sum is empty).

### 5.5 Strong-coupling `SU(N)` 4D — the one positive anchor, and it is P3 `[T]` (mod `[R]`)

At `β<1/48`, the dual-surface expansion (1.1) is absolutely convergent, its minimal tube is the unique
leading class, and the recoupling positivity gives `f₁²≥r₀(β)>0`. This is *exactly* P3 in surface
language (`P3-FOURWALL`), and it is the **only** rigorous positive binding bound in the vortex picture.
It is **not** a vortex-percolation result (the positivity is `1/N>0`, §4.1), and it gives **no**
isolation `ρ` and does not transport past the `1/48` Ursell wall.

| test | confining? | vortex picture rigorous? | `f₁²(tr F²)` | `ρ` | lesson |
|---|---|---|---|---|---|
| 4D `Z(N)` | confined phase | **exact** (positive dual, FK) | `>0` | computable, not closed form | abelian luxury |
| 3D `Z(N)`/`U(1)` | yes (monopole) | **exact**, semiclassical | `>0` | computable | abelian+3D |
| 2D YM | **yes** (`σ>0`) | exact (Migdal) | **`0`** | `0`/degenerate | **`σ>0 ⇏` atom** |
| free field | no | trivial | `0` | `1` (none) | acid pass |
| SU(N) 4D `β<1/48` | yes (SZZ) | convergent dual = P3 | `≥r₀>0` | not given | binding = recoupling |
| SU(N) 4D `β≈2.5` | `[O]` | **signed → none** | `[O]` | `[O]` | **the wall** |

---

## 6 · What the vortex picture gives, and the exact missing estimate `[O]`

### 6.1 What it gives

- **A rigorous confinement inequality in the charged sector** `[T]`: Tomboulis–Yaffe `⇒ σ>0` (given
  vortex condensation). At β<1/48 (SZZ) / β<12/17 (Borgs–Seiler) the premise holds; at β≈2.5 it is `[O]`.
- **An exact dual representation of the `tr F²` correlator** `[T]` (1.1), whose strong-coupling
  truncation *is* P3 (`f₁²≥r₀>0`, `[C]` mod `[R]`).
- **A proof that the central vortex object is neutral-trivial** `[T]` (Lemma CV-1), hence that the
  vortex *condensate* cannot carry the atom.

### 6.2 The exact missing estimate `(E-CV)` `[O]`

The bound the vortex picture *would* need to give a binding result is a positive, volume-uniform
lower bound on the neutral-sector dual-surface amplitude, with sign control over the remainder:

> **`(E-CV)`.** There exist `c_v>0`, and `a₀,β₁,L₀` such that, uniformly for `β` through `≈2.5`,
> `a<a₀`, `L>L₀`, the connected `tr F²` surface sum obeys
> ```
>   | Σ_{S: ∂S=P_x∪P₀} W_v(S) |  ≥  c_v ( Σ_{S} |W_v(S)| ) ,       W_v = Weingarten-signed,      (E-CV)
> ```
> i.e. a **sign-cancelling positive regrouping** of the dual-surface sum in the **N-ality-0** sector
> into a Gibbs/FKG measure whose connection probability lower-bounds the tube amplitude — equivalently,
> the **dressed-tube resummation lower bound** of `P3-FOURWALL §3.2 [R]` extended past `β=1/48` **with
> sign control**, uniformly in `(L,a)`.

Three honest namings, and they are one object:
- **vortex:** a positive (FK/Edwards–Sokal) measure for the **neutral** surface ensemble at β≈2.5 —
  absent because the Weingarten weights are signed (`PROOF-RANDOM-CLUSTER`, `[NG]`), and the centre
  phase that would enforce positivity is **1** on the neutral sector (Lemma CV-1).
- **binding:** the positive lower bound `f₁²≥r₀>0` with isolation `ρ≥1+c` — R2's `(E-I2)`, the order-2
  `TP₂` datum (`[O]`).
- **random surface:** a sign-cancelling regrouping of the `SU(N≥2),d=4` surface sum (the open sign
  problem flagged by Cao–Nissim–Sheffield, `PROOF-CONFINEMENT-RANDOM-SURFACE §2–3`).

**They coincide.** That is the finding: the centre-vortex bound on binding is *literally* the same
missing estimate as `(E-I2)` and as the non-abelian random-surface sign problem — not a new, weaker,
or more tractable object. The vortex representation renames the wall; it does not lower it.

### 6.3 The conditional (why granting it is circular) `[C]`

> **`(CV-★)` `[C]`.** Grant (i) `(E-CV)` (neutral positive surface measure at 2.5), (ii) a vortex–glueball
> vertex bound `|⟨0|tr F²|G⟩|² ≥ (c_v ρ̂)²`, and (iii) `m_{0⁺⁺}=c_g√σ` with `c_g` bounded. Then
> `f₁²≥r₀>0` and `ρ≥1+c` follow. **But** (ii) is the neutral binding content and (iii) is the gap—string
> transmutation; the conditional re-states the target rather than deriving it. No content is added; the
> vortex picture buys nothing here over `(T★)`.

---

## 7 · Status ledger

| item | statement | label |
|---|---|---|
| centre-vortex phase on N-ality-0 `R` | `χ_R(zI)/d_R = ω^{q}=1` for `q=0`; `=ω^q` else | **[T]**+**[V]** `[A]` |
| `tr F²` neutral-blindness | `tr F²` is N-ality 0 ⇒ every vortex phase is 1 | **[T]** |
| dual expansion of `⟨tr F² tr F²⟩_c` | `Σ_S W_v(S)`, `W_v` Weingarten-signed = CPS rep. | **[T]** (representation) |
| minimal tube in dual picture | minimal spanning cylinder `4N_t`, coincides with 4-wall tube; N-ality 0 | **[T]** (`P3-FOURWALL`) |
| center projection neutral sector | identically trivial (`max|W_A−1|=0`; adjoint `=1`) | **[T]**+**[V]** (`probe-center`) |
| Tomboulis–Yaffe | `⟨W_fund⟩≤2[½(1−Z⁻/Z)]^{…}` ⇒ `σ>0` (charged only) | **[T]** |
| `σ>0 ⇏ tr F²` atom | 2D YM: `σ_fund>0` yet `f₁²(tr F²)=0` exactly | **[T]**+**[V]** `[B]` |
| `m∝√σ` bridge | measured `c≈3.5–4`, not a theorem; area law ⇏ `0⁺⁺` gap | **[O]** (empirical) |
| strong-coupling tube positivity | `a₀>0` from recoupling `1/N>0`, `u>0` (P3), NOT percolation | **[T]** `[C]` mod `[R]` |
| P3 residual gap | dressed resummation `[R]`: lower bound from signed remainder | **[R]** |
| abelian vortex picture exact | `Z(N)`4D/3D: positive weights, FK/FKG gap ⟺ percolation | **[T-standard]**+**[V]** |
| non-abelian neutral measure | Weingarten-signed (Euler char SU(2); ~half<0 SU(N)) ⇒ no FKG | **[NG]** |
| overlap lower bound `(4.2)` | collapses at 2.5 (signed `W_v`); `≥0` only | **[NG]** |
| isolation `ρ` from vortices | needs order-2 second cone on closed surfaces — absent | **[NG]**/`[O]` |
| `(E-CV)` = `(E-I2)` = random-surface sign | sign-cancelling positive regrouping at 2.5 | **[O]** |
| `(CV-★)` conditional | granting it is circular (restates target) | **[C]** |
| **Clay mass gap** | unbuilt | **[O]** · **0%** |

---

## 8 · Witness — `research/v3/w7/t5_center_vortex.py`

`py -3.12 -I t5_center_vortex.py` (exit 0; numpy 2.4.6; deterministic quadrature, no MC, ~1 s).
Quoted output (abridged):

```
[A] SU(2) -I: fund χ/d = -1.000000 ; adjoint χ/d = +1.000000 ; j=2 χ/d = +1.000000
    SU(3) k=1: fund phase = -0.500000+0.866025j ; adjoint(N-ality 0) phase = +1.000000
    max |neutral phase − 1| = 4.4e-16                          => tr F^2 is centre-blind  [T]
[B] Haar orthogonality ∫χ_j dU = 1.000, 0.000, 2.9e-17, 2.9e-17, -7.1e-18   (δ_{j0})
    β=2.5: σ_fund = 0.60000 > 0 (CONFINES) ; f1^2(tr F^2) = 0  (same theory)     [T]+[V]
[C] u(β)=I2/I1 : 3.50e-3 (β=.01), 5.69e-3 (β=1/48), 0.1237 (β=.5) ; recoupling 1/N=1/2>0
[D] free field: 2-particle cut, point-mass weight f1^2 = 0 (ACID PASS)
[E] fundamental phase drives σ (charged); adjoint/tr F^2 phase = +1 (neutral-blind); Weingarten signed [NG]
```

All quoted numbers are reproduced exactly by the script; the group-theoretic identities
(Lemma CV-1), the 2D orthogonality (`[B]`), and the free-field cut (`[D]`) are `[T]`, the witness is a
self-check `[V]`, and the Weingarten-sign statement is `[T-standard]` cited from
`PROOF-RANDOM-CLUSTER`. No `[C]` upgraded to `[T]`.

---

## 9 · Sharpest line

**The centre-vortex picture bounds confinement, not binding.** Its one rigorous inequality
(Tomboulis–Yaffe `⇒ σ>0`) lives in the **charged** N-ality-1 sector; `tr F²` is a colour singlet, and
the centre phase it sees is `ω⁰=1` **exactly** (Lemma CV-1, `[A]`), so the entire vortex-condensate
mechanism is provably neutral-blind. That `σ>0` does not give the atom is a theorem-by-example, not a
gap: 2D pure YM confines (`σ_fund>0`) with `f₁²(tr F²)=0` exactly (`[B]`), and the `m_{0⁺⁺}∝√σ` bridge
is measured, not proved. The only route that could bound the tube overlap — a **positive** dual-surface
(vortex) measure — exists for `Z(N)/U(1)` (FK/FKG) and is **sign-killed** for non-abelian `SU(N)` in
4D (Weingarten; Euler-characteristic signs). So the attempted overlap bound `(4.2)` collapses to `≥0`
at β≈2.5, and the isolation `ρ` needs the same order-2 second-cone that the gauge quotient lacks. The
**exact missing estimate** is a sign-cancelling positive regrouping of the dual-surface weight in the
N-ality-0 sector at β≈2.5 — which is **literally** `(E-I2)` and the non-abelian random-surface sign
problem, not a weaker object. **The vortex representation renames the wall; it does not lower it. No
QED. Clay 0%.**

---

*Refs.* **Corpus (read):** `R-WAVE-SYNTHESIS.md`, `R5-EXACT-ANCHORS.md` (`§1.1` 2D `f₁²=0`; `§1.4`
P3; `§3` transport), `R2-ISOLATION.md` (`§3.3` order-2; `(E-I2)`), `docs/paper-v2/ch04-gap.tex`;
`unconditional/PROOF-GAUSS-INVARIANT-PF.md`, `PROOF-CENTER-VORTEX.md` (TY `σ>0`, `ρ̂` at threshold),
`PROOF-RANDOM-CLUSTER.md` (Weingarten signs, FK abelian-only), `PROOF-CONFINEMENT-RANDOM-SURFACE.md`
(CPS/CNS small-λ; sign problem), `P3-FOURWALL-DERIVATION.md` (four-wall tube, `a₀>0`, `[R]`
resummation); `v3/P3-POLYAKOV-CENTER.md` (centre-blind constant mode); `research/probe-center/REPORT.md`
(`§3–4` neutral-blind, `W_A(proj)=1`). **Standard inputs (textbook, not re-fetched):** Tomboulis–Yaffe,
CMP **100** (1985) 313; 't Hooft, NPB **138** (1978) 1, NPB **153** (1979) 141; Brydges–Giffen–
Durhuus–Fröhlich, NPB **275** (1986) 459; Cao–Park–Sheffield, arXiv:2307.06790 (Comm. AMS);
Cao–Nissim–Sheffield, arXiv:2505.16585, arXiv:2509.04688; Fortuin–Kasteleyn, Physica **57** (1972) 536;
Edwards–Sokal, J. Stat. Phys. **52** (1988) 223; Wegner, J. Math. Phys. **12** (1971) 2259; Migdal,
Sov. Phys. JETP **42** (1976) 413; Morningstar–Peardon, PRD **60** (1999) 034509. **Witness:**
`research/v3/w7/t5_center_vortex.py` (run once, `py -3.12 -I`, numpy only). **No third-party sources
downloaded; no `[C]` upgraded to `[T]`; no QED; Clay 0%.**
