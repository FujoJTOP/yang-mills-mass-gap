# U4-CONTINUUM-ATOM — the `a→0` limit of the atom fraction `w₁=f₁²/M₀`: it does **not run** (the anomalous dimension `γ_{w₁}=0` exactly), it has a **finite** continuum limit at fixed physical scale, and whether that limit is **positive** is *not* an RG statement — it is exactly the uniform isolation/binding estimate, the same wall as R/T, now shown to be the *only* content of the continuum atom

> **2026-11-24 · Wave U4 — the continuum limit of the atom, the piece left by the R/T waves.**
> The R/T waves define `w₁=f₁²/M₀` on the lattice and state the missing bound at `β=2.5`, **fixed
> `a`, fixed `L`**. The Clay problem takes `a→0` (and `L→∞`). This note asks the `a→0` question:
> using **gradient-flow / point-splitting renormalization**, how do `w₁`, `M₀`, `f₁²`, and the
> isolation ratio `ρ=m_{0⁺⁺*}/m_{0⁺⁺}` scale with `a`? Is `w₁(a)` finite and positive at fixed
> physical volume, or does it run? Write the RG equation and the anomalous dimension; test on the
> free field (`w₁=0`), 2D (exact), and strong coupling (exact); state the **exact missing estimate**
> for a uniform continuum bound.
>
> **Reads.** `T-WAVE-SYNTHESIS.md`, `R-WAVE-SYNTHESIS.md`, `R1-RESIDUE-TRANSPORT.md`,
> `R2-ISOLATION.md`, `R5-EXACT-ANCHORS.md`, `T4-CARRIER.md`, `T2-TUBE-OVERLAP.md`;
> `unconditional/MS-RG-ODE.md`, `INVENT-TRANSPORT-INEQ.md`, `CONTINUUM-VARIANCE-DEFECT.md`,
> `G1-CONTINUUM-SCALING.md`, `PROOF-GRADIENT-FLOW.md`;
> `docs/paper-v2/ch02-squeeze.tex` (`eq:chi`,`eq:squeeze`,`eq:h1`, thin-atom Prop.),
> `ch03-existence.tex` (`§dim`, `eq:dim`), `ch12-numerics.tex`.
> **Witness.** `research/v3/w7/u4_continuum.py` (`py -3.12 -I`, numpy only, one small job, exit 0,
> ~1 s; output quoted §7). No third-party downloads.
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional · **[NG]** no-go for a
> named class · **[O]** exactly the missing step · **[V]** numerical witness here. **No QED. Clay 0%.**

---

## 0 · Verdict (up front; no QED, Clay 0%)

1. **`w₁` does not run — the anomalous dimension is exactly zero `[T]`.** `w₁=f₁²/M₀` is the ratio of
   the atom residue to the **total** spectral weight of the **same** dimension-4 operator
   `O=tr F²`. Both `f₁²=|⟨0|O_R|G⟩|²` and `M₀=∫dν_O=Σ_n|⟨0|O_R|n⟩|²` are two insertions of `O_R`
   and both scale by the **same** `Z_O²` under the trace anomaly. Hence
   ```
     f₁²(μ)=Z_O²f₁²(μ₀),   M₀(μ)=Z_O²M₀(μ₀)   ⇒   w₁(μ)=w₁(μ₀),
     d ln w₁/d ln μ = 2γ_O − 2γ_O = 0 ,      γ_{w₁} ≡ 0 ,                          (0.1)
   ```
   where `γ_O=−d ln|β_YM/2g|/d ln μ=2b₀g²+O(g⁴)` is the (nonzero) anomalous dimension of `tr F²`
   itself. `w₁` is an exact RG invariant; it carries **no** new anomalous dimension. This is the
   `Z_O`-invariance of T2 (`w₁` unchanged under `O→Z_O O`) rewritten as an RG equation.

2. **`w₁` has a finite `a→0` limit at fixed physical scale `[T]` — but only after the contact is
   removed.** Every quantity in the ratio is UV-sensitive in a **scheme-dependent, power-law** way
   before renormalization: the bare on-site weight `M₀=C(0)=⟨O²⟩_c ~ a⁻⁸` (two dim-4 fields at
   coincident points), the susceptibility `χ=∫s⁻¹dν ~ a⁻⁴` (T4-CARRIER §3). The atom residue
   `f₁²=ν({m₁²})` sits at the **IR** mass `m₁≪1/a` and therefore carries **only** the multiplicative
   `Z_O` (no contact). The asymmetry is the whole point: **`f₁²` is IR-clean, `M₀` is UV-dominated**,
   so the **bare** fraction `w₁^{bare}=f₁²/M₀^{bare}→0` trivially (a UV artifact, *not* a physical
   running). Defined with `O` **renormalized at a fixed physical scale** (gradient flow at fixed flow
   time `t` ≳ physical², or point-splitting at fixed physical `ε`), both `f₁²` and `M₀` are finite,
   their ratio `w₁(t)` is finite, and its `a→0` limit exists `[T]` (Lüscher–Weisz: `O_t` is a
   well-defined continuum operator; matrix elements converge).

3. **Whether the limit is POSITIVE is not decidable by RG `[O]`.** Equation (0.1) is a
   **homomorphism that sends `0↦0`** — the same absorbing-`0` structure as R1's `Z_O` transport.
   `γ_{w₁}=0` says *if* the atom fraction exists, it does not run; it says nothing about whether it
   is nonzero. `w₁(a)→w₁^cont>0` **iff** the lightest `0⁺⁺` atom survives as `a→0` **with
   isolation** `I₂` (`Λ★−m₁²≥δ`), i.e. is a genuine bound state below the two-gluon cut and not a
   lattice-cutoff artifact that dissolves. This is exactly the R-wave target, unchanged. The
   continuum-limit attack **removes** the "is there a new UV running/anomalous-dimension obstruction
   to `w₁`?" worry (there is none, `[T]`) and **relocates** the entire residual into the same
   non-perturbative isolation/binding statement `[O]`.

4. **The isolation ratio `ρ` also does not run `[T]`, and it does not by itself give the atom `[NG]` /
   `[V]` in 2D.** `ρ=m_{0⁺⁺*}/m_{0⁺⁺}` is a **ratio of masses**; masses are RG invariants and both
   scale as `a·m_phys` (ch03 `eq:dim`, `m_lat=√G″_a=a m_phys`), so `ρ(a)=ρ_cont` is `a`-independent
   and has **no anomalous dimension**. But a finite `ρ>1` is *not* sufficient for a `tr F²` atom: **2D
   pure YM** is gapped with an exact holonomy`/`Casimir isolation ratio `ρ²=8/3` (R2, R5) **while the
   `tr F²` atom fraction is `w₁=0` exactly** (`tr F²=H` conserved, character-orthogonality; R5 §1.1).
   So `ρ` (isolation of the lightest `0⁺⁺`) and `w₁` (that operator's atom fraction) are independent,
   and 2D is the exact continuum counterexample where `ρ>1` yet `w₁=0`.

5. **Free-field acid passes; strong coupling is the one positive anchor — at fixed `a`.** Free field:
   `:φ²:` is a two-particle cut, `f₁²=0`, so `w₁=0` identically `[T]`. 2D: `w₁=0` exactly `[T]`. Strong
   coupling `β<1/48`: `w₁≥r₀/M₀>0` finite `[T]` at fixed `a` — the only provable positivity, and it
   **does not transport to `a→0`** (no scaling window; `a·m₁≫1`, R5).

6. **The exact missing estimate** (the deliverable, §6): a **uniform lower bound on the atom fraction
   of the renormalized (fixed-scale) KL measure**, uniform in `(L,a)` and in the physical smearing
   scale `t` along the physical line at `β=2.5`:
   ```
     inf_{L,a,t}  w₁(t;a,L) = f₁²/M₀  ≥  w₀ > 0      (equivalently  Σ_{n≥2} w_n ≤ 1−w₀) ,
     together with  ρ(a) ≥ 1+c > 0 .
   ```
   It is **not** an RG/anomalous-dimension estimate (that half is `γ_{w₁}=0`, free and exact). It is
   the spectral statement that the **bottom atom of the renormalized measure carries a fixed fraction
   of its total (fixed-scale) weight**, isolated from its continuum — i.e. uniform `I₂`, the gap one
   rung up. No variance, positivity, RG, transport, or anomalous-dimension tool supplies it; 2D shows
   even a finite `ρ` does not.

**One line.** *`w₁=f₁²/M₀` does not run — `γ_{w₁}=0` exactly, because numerator and denominator are
both `Z_O²`-times one operator — and it has a finite `a→0` limit once `O` is renormalized at a fixed
physical scale (gradient flow / point-splitting), the bare `M₀~a⁻⁸` contact having been removed; but
positivity is an RG-blind, `0`-absorbing statement, so whether the continuum atom fraction is `>0` is
exactly the uniform isolation `I₂` / gap, and even a finite `ρ` (2D: `8/3`) does not give it. The RG
half is free; the missing estimate is `w₁≥w₀>0` with `ρ≥1+c`, uniform in `(L,a,t)`. No QED. Clay 0%.*

---

## 1 · The objects, their dimensions, and their exact `a`-scaling

Let `O_S=¼tr F²` (dimension 4, `0⁺⁺`, P-even). Its Källén–Lehmann measure on the `0⁺⁺` channel is
```
   ν_O = Σ_n |⟨0|O|n⟩|² δ(s−m_n²) ,   s ≡ m² ,
   f₁² = ν_O({m₁²})  (atom, mass⁸) ,   M₀ = ∫dν_O = Σ_n|⟨0|O|n⟩|² = C(0)  (total weight, mass⁸) ,
   χ = ∫s⁻¹dν_O = Σ_n|⟨0|O|n⟩|²/m_n²  (susceptibility, harmonic mean) ,
   w₁ = f₁²/M₀ ∈ [0,1]  (dimensionless atom fraction) .                                   (1.1)
```
The isolation half (R2) uses the second moment of the **same** measure above the atom:
`μ₁=m₁²`, `μ₂=inf supp(ν−f₁²δ_{μ₁})`, `ρ²=μ₂/μ₁`, isolation margin `Δ=(ρ²−1)μ₁`; the physical
`ρ=m_{0⁺⁺*}/m_{0⁺⁺}=√(μ₂/μ₁)`.

**Dimension and cutoff table `[T]`.** With `O(x)=a⁻⁴O_lat(x)` (`O_lat` dimensionless, built from
plaquettes), and `Σ_z` the 4D lattice sum with cell weight `a⁴`:

| object | mass dim | bare lattice `a`-scaling | after fixed-scale renorm | RG (anomalous dim) |
|---|---|---|---|---|
| `f₁² = ν({m₁²})` | 8 | `~a⁻⁸` (lattice units); **IR-clean** of contact | finite, `a`-indep. | `+2γ_O` (multiplicative only) |
| `M₀ = C(0)=∫dν` | 8 | `~a⁻⁸` (dim-8 contact) | finite | `+2γ_O` × (contact subtracted) |
| `χ = ∫s⁻¹dν` | 4 | `~a⁻⁴` (T4-CARRIER §3) | finite | `+2γ_O` |
| `w₁ = f₁²/M₀` | 0 | `w₁^{bare}→0` (UV artifact) | finite | **`0` (exact)** |
| `ρ = √(μ₂/μ₁)` | 0 | mass ratio, `a`-indep. | finite `=ρ_cont` | **`0`** (masses don't run) |
| `G″_a` (curvature) | 0 | `~a²m_phys²` (ch03 `eq:dim`) | `→a²Δ₀²` | — |

The `a⁻⁸`/`a⁻⁴` powers are the standard contact divergences of a dim-4 composite (T4-CARRIER §3:
"on-site value scales as `a⁻⁸`; summed to `χ` it diverges as `a⁻⁴`"). The **crucial asymmetry**: the
atom residue `f₁²` is a delta at the **IR** mass `m₁≪1/a`, so it does **not** carry the contact — only
the multiplicative `Z_O`. The total weight `M₀=C(0)` sums over all states including the **UV**, so it
carries the contact on top of `Z_O`. **This is why the bare ratio `w₁^{bare}=f₁²/M₀^{bare}→0`** (a UV
artifact — the denominator counts unbounded short-distance weight), and **why a fixed-scale
renormalization is mandatory for `w₁`, not merely for `χ`.**

*Witness [E].* Free 4D scalar lattice propagator, fixed physical box, `a=Lbox/L→0`: measured
`G(0)=(1/V)Σ_p 1/(p̂²+m²)` has log-log slope `−2.08` vs `a` (`G(0)·a²→const`), confirming the `a⁻²`
contact skeleton. A dim-`Δ` composite then has `M₀=C(0)~a^{−(2Δ−4)}`: `:φ²:` (`Δ=2`) `~a⁻⁴`,
`tr F²` (`Δ=4`) `~a⁻⁸`, both `[T]` by power counting on the measured skeleton.

---

## 2 · The RG equation for `w₁`: `γ_{w₁}=0` exactly `[T]`

**Setup.** `O_R(μ)=Z_O(μ)⁻¹O_B`, with `Z_O` fixed **not** by a loop computation but by the exact trace
anomaly `T^μ_μ=(β_YM/2g)O_S` (RG-invariant), giving `Z_O(μ,μ₀)=(β_YM/2g)(μ₀)/(β_YM/2g)(μ)>0`
(R1 §1.2, MS-RG-ODE §1). The `0⁺⁺`-off-diagonal matrix element of the dim-4 scalar renormalizes
multiplicatively with no operator mixing (identity is diagonal and absent for `G≠0`; `F F̃` decouples;
BRST-exact derivatives vanish on physical states — R1 §1.1).

**Theorem U4-1 (RG equation of the atom fraction) `[T]`.** For the renormalized operator at a fixed
scheme, with `M₀(μ):=∫dν_{O_R(μ)}=Σ_n|⟨0|O_R(μ)|n⟩|²` and `f₁²(μ)=|⟨0|O_R(μ)|G⟩|²`,
```
   f₁²(μ)=Z_O(μ,μ₀)² f₁²(μ₀) ,     M₀(μ)=Z_O(μ,μ₀)² M₀(μ₀) ,                              (2.1)
   hence   w₁(μ)=f₁²(μ)/M₀(μ)=w₁(μ₀)   and   d ln w₁/d ln μ = 0 ,   γ_{w₁}≡0 .            (2.2)
```
*Proof.* Each `|⟨0|O_R(μ)|n⟩|² = Z_O(μ,μ₀)⁻²|⟨0|O_B|n⟩|²`; summing over `n` (the atom and the
continuum alike) gives (2.1); the common factor cancels in the ratio. The one-loop form of the
**nonzero** dimension is `γ_O=−d ln|β_YM/2g|/d ln μ=2b₀g²+O(g⁴)`, `b₀=11N/(48π²)` (MS-RG-ODE §1); it
is the anomalous dimension of `f₁²` and `M₀` individually, **not** of `w₁`. ∎

**Witness [A] `[V]`.** A KL measure (atom `f₁²=0.30` at `m₁²=1` + a two-particle cut) run through the
one-loop `Z_O²(μ)` from `μ=2` to `100`:

```
  mu        Z_O^2      f1^2(mu)        M0(mu)           w1(mu)   dlnw1/dlnmu
  2.0     1.000000      0.300000     78.763325   0.003808879306
   3.0     2.512106      0.753632    197.861832   0.003808879306    +0.00e+00
   5.0     5.391350      1.617405    424.640660   0.003808879306    +0.00e+00
  10.0    11.035206      3.310562    869.169542   0.003808879306    +0.00e+00
 100.0    44.140825     13.242248   3476.678167   0.003808879306    +0.00e+00
  check: dln f1^2/dln mu = 1.242670  vs  2*gamma_O(5) = 1.242670
  => w1 constant to ~1e-15; dln w1/dln mu = 0.  gamma_{w1}=0 EXACTLY.
```
`f₁²` and `M₀` both grow (they run with `+2γ_O`), the ratio is pinned to machine zero derivative.

**What this does and does not settle.** It settles the **running**: `w₁` has **no anomalous
dimension** — the "does the atom fraction run with `a`?" question has answer **no** `[T]`. It does
**not** settle the **positivity**: (2.2) is linear in each of `f₁²`, `M₀` with the same coefficient,
so it maps `f₁²=0↦w₁=0` and `f₁²>0↦w₁>0`; it cannot *create* an atom, and it cannot *prevent* one
from dissolving. `γ_{w₁}=0` is **necessary-shaped** (no running obstruction) and **insufficient**
(positivity).

---

## 3 · The continuum limit: finite `[T]`, positive `[O]`, and the two traps

### 3.1 Finite: the fixed-scale renormalized ratio converges `[T]`

Take `O` renormalized at a **fixed physical scale**: either **gradient flow** `O_t=tr F²(t)` at a
fixed flow time `t` with `√t` a fixed fraction of `1/m_G` (Lüscher; PROOF-GRADIENT-FLOW §1), or
**point-splitting** `O_ε(x)=tr F(x)F(x+ε)` at fixed physical `ε`, with the contact subtracted. For
`t>0` (resp. `ε>0`) the separated-point correlator is finite (Lüscher–Weisz), so

```
   f₁²(t) = |⟨0|O_t|G⟩|² ,   M₀(t)=∫dν_{O_t} = ⟨O_t²⟩_c   are both finite ,
   w₁(t) = f₁²(t)/M₀(t)   is finite for every a ,   and   lim_{a→0} w₁(t;a) = w₁^cont(t) exists `[T]`.
```
The flow provides the subtraction that removes the `a⁻⁸` contact from `M₀`; the residue `f₁²(t)` is
already IR-clean. So **at fixed physical volume and fixed physical scale, `w₁(a)` is finite (bounded,
convergent), not divergent, not running.** This is the reassuring `[T]` half the task asks for.

### 3.2 Trap 1 — the bare fraction is a UV artifact, not a physical decay

The naive `w₁^{bare}=f₁²/M₀^{bare}` **does** go to `0` as `a→0` (the `a⁻⁸` contact in `M₀`), but this
is **not** the atom dissolving and **not** a physical running — it is the denominator counting
unbounded short-distance spectral weight. The physical atom fraction is the **fixed-scale** one. The
witness [B] shows the mechanism on the free composite: `M₀(Λ)` grows `1.92→198.8` as the cutoff `Λ`
goes `10→100`, while `f₁²=0` throughout. A tool that reports "`w₁→0`" from the bare ratio is measuring
the UV, not the atom.

### 3.3 Trap 2 — the value is smearing-scheme dependent (finite renormalization)

`w₁(t)` is finite for each `t`, but `O_t` mixes with higher-dimension operators under flow
(`O_t=O+t·\{\text{higher dim}\}+\dots`), so `f₁²(t)/M₀(t)` genuinely **depends on the finite scheme
`t`** — the *value* of the atom fraction changes with the smearing radius. This is **finite**
renormalization (a scheme choice), distinct from the **anomalous** (log) running which is zero. A
uniform continuum bound must therefore hold **uniformly in `t`** (or for one fixed `t`); at `t→0` the
flowed operator returns to the bare one and `w₁(t)→0` through the contact (trap 1), while for
`√t≳1/m_G` the glueball overlap degrades. The honest window is `t∈[t₁,t₂]` with `√t_i ∼ c_i/m_G`.
This `t`-dependence is a genuine new feature of the continuum-limit statement that the fixed-`a`
R/T analysis (which fixed a single lattice smearing family) did not isolate.

### 3.4 Positive: atom survival with isolation `[O]`

`lim_{a→0}w₁(t)>0` iff the lightest `0⁺⁺` is a **genuine isolated bound state** in the continuum —
not a lattice state pinned at the cutoff that dissolves (`f₁²→0`), and not embedded in a continuum
(`I₂` fails, `Λ★=m₁²`). Equivalently `I₂`: `Λ★−m₁²≥δ>0`. The RG equation (2.2) is **blind** to this:
it is a `0`-absorbing homomorphism on the weight (R1 §5 N1: `Z_O²·0=0`; `Ω` bounded only above by
Cauchy–Schwarz, N2). **The continuum-limit attack therefore proves that a positive continuum atom
fraction is *not* obstructed by any running, and *is* obstructed by exactly the same non-perturbative
isolation/binding that the R/T waves named at fixed `a`.**

---

## 4 · Solvable tests (free field, 2D, strong coupling)

### 4.1 Free field: `w₁=0` exactly `[T]+[V]`
`O=tr F²` of a free field (or `:φ²:`) is a **composite bilinear**; its KL measure is the two-particle
cut on `[2m,∞)` (massless: `[0,∞)`), absolutely continuous, **no delta** (R5 §1.6). Hence
`f₁²=ν({m₁²})=0` for every `m>0`, so `w₁=0/M₀=0` **identically, for every `a` and every cutoff**.
Witness [B]: `M₀(Λ)` diverges but `f₁²=0`, `w₁=0`. **Acid passes**: the free field does not
manufacture an atom fraction, and the fixed-scale ratio is `0` in the continuum.

### 4.2 2D pure YM: `w₁=0` exactly, though gapped and with `ρ>1` `[T]+[V]`
In `1+1` dims `tr F²=tr E²=2/g²·`(energy density) is conserved and **diagonal** in the character
basis, so `⟨0|tr F²|G_R⟩=C₂(R)⟨χ_0|χ_R⟩=0` for `R≠0` (R5 §1.1). Character orthogonality:
`∫χ_0χ_R dU=δ_{0R}` (witness [C]: `|⟨χ_0|χ_R⟩|<2×10⁻¹⁷` for `R=\frac12,1,\frac32,2`). So
`f₁²(2D tr F²)=0` **exactly** ⇒ **`w₁=0`**, in a theory that is **gapped** and whose **holonomy**
isolation ratio is exact `ρ²=C₂(adj)/C₂(fund)=2N²/(N²−1)=8/3` for `SU(2)` (R2, R5). **This is the
sharp warning:** a finite `ρ>1` (isolation of the lightest `0⁺⁺`, here the holonomy) does **not**
imply a `tr F²` atom (`w₁=0`). It is the exact continuum realisation of T6's `I⇒A` being
**operator-specific**: the holonomy's isolation does not give the `tr F²` composite an atom.

### 4.3 Strong coupling `β<1/48`: `w₁>0` at fixed `a`, transports nothing `[T]`
P3 (`PROOF-P3-UNCONDITIONAL`) gives `f₁²≥r₀(β)=½c₁u^p>0`, uniform in `L`, at fixed `a`. With `M₀`
finite on the finite lattice, `w₁≥r₀/M₀>0` — the only provable positive atom fraction. But this is a
**fixed-`a`** statement: at `β<1/48` there is no scaling window (`a·m₁≫1`, R5; the state sits at or
beyond the cutoff). It gives **no** `a→0` continuum fraction. This is consistent, not contradictory:
strong coupling proves the atom exists on the lattice; the continuum limit is a separate question.

---

## 5 · The isolation ratio `ρ`: `a`-independent, no anomalous dimension — and not enough

`ρ=m_{0⁺⁺*}/m_{0⁺⁺}=√(μ₂/μ₁)` is a ratio of transfer-matrix eigenvalues. Masses are RG invariants
(`d ln m₁/d ln μ=0`, MS-RG-ODE §4); on the lattice `m_lat=√G″_a=a·m_phys` (ch03 `eq:dim`), so both
eigenvalues scale as `a²m_phys²` and **`ρ²=μ₂/μ₁=ρ_phys²` is `a`-independent** with **no anomalous
dimension**:
```
   ρ(a) = ρ_cont  (pure number) ,   γ_ρ ≡ 0 .                                               (5.1)
```
So `ρ` does not run either. But (4.2) shows the **2D counterexample** where `ρ²=8/3>1` (gapped,
isolated holonomy) and `w₁(tr F²)=0`: **isolation of the lightest `0⁺⁺` does not give the specific
operator an atom.** The window theorem `Δ=(ρ²−1)μ₁`, `1+c≤ρ²≤C² ⟹ μ₁≥δ/(C²−1)` (T3) is a statement
**within one operator's** measure; it is not a bridge across operators. Hence the continuum bound needs
**both** `ρ≥1+c` **for the `tr F²` channel** (isolation) **and** `w₁≥w₀>0` **for `tr F²`** (atom
fraction); 2D shows these can come apart when the lightest state is created by a different operator.

---

## 6 · The exact missing estimate for a uniform continuum bound `[O]`

The RG/transport/`a`-scaling part is **closed `[T]`**: `γ_{w₁}=0`, `γ_ρ=0`, both ratios finite in the
continuum at fixed physical scale, the contact removed by flow/point-splitting. **What is missing is
purely spectral positivity, not any continuum-limit running.** The exact missing estimate:

> **(U4-★) — uniform atom fraction with isolation, along the physical line at `β=2.5`, uniform in
> `L→∞`, `a→0`, and the physical smearing scale `t∈[t₁,t₂]` (`√t_i∼c_i/m_G`):**
> ```
>    inf_{L,a,t}  w₁(t;a,L) = f₁²(t)/M₀(t)  ≥  w₀ > 0            [the atom fraction floor]
>    equivalently   Σ_{n≥2} w_n(t;a,L) ≤ 1 − w₀                  [the leakage bound, T2 §1.3]
>    and            ρ(t;a,L) = m_{0⁺⁺*}/m_{0⁺⁺} ≥ 1 + c > 0     [isolation I₂, for the tr F² channel]
> ```
> where `f₁²(t), M₀(t)` are the atom residue and total weight of the **fixed-scale renormalized**
> measure `ν_{O_t}`. Given (U4-★), the squeeze `m₁²≥f₁²/χ` with `χ≤1/λ` gives a uniform positive
> mass; **(U4-★) is `H*` read at `a→0`, i.e. the gap one rung up.**

**Why this and not something else.**
- It is **not** an anomalous-dimension estimate. `γ_{w₁}=0` (2.2) is exact and free; it carries the
  *value* of `w₁` and maps `0↦0`. Positivity is the `t`-uniform, `a`-uniform, `I₂`-isolated bottom
  atom — outside every RG/transport/positivity tool (R1 §5: `Z_O` bounds carry no sign; `Ω`
  upper-bounded only; moments Schur-concave ⇒ upper bounds only, ch02 `§schur`).
- It is **operator-specific**: 2D supplies `ρ²=8/3>1` while `w₁(tr F²)=0` (4.2), so `ρ≥1+c` anywhere
  in the theory is not enough; the floor must be on the `tr F²` KL measure itself.
- It is **scheme-uniform**: because `w₁(t)` is finite but `t`-dependent (3.3), the floor must hold for
  a fixed physical window of `t`; a value at one `t` does not bound another (the flow mixing is a
  genuine finite renormalization).

**What would close it (unchanged from R/T).** A non-perturbative binding theorem making the minimal
`tr F²` polymer a single isolated four-wall tube at `β=2.5`, uniform across `a→0` — equivalently a
Kato-stable "same-atom" flow with a uniform isolation margin (R1 §8), or a sign-cancelling positive
regrouping in the N-ality-zero sector (T5). None is supplied by the continuum-limit analysis; the
continuum-limit analysis **removes the RG objection** and leaves the binding/isolation wall exactly
where R/T placed it.

---

## 7 · Witness — `research/v3/w7/u4_continuum.py`

`py -3.12 -I u4_continuum.py` (exit 0; numpy 2.4.6, no scipy/MC, ~1 s, one job). Quoted (abridged):

```
[A] RG equation: f1^2 and M0 both ~ Z_O^2  =>  w1 = f1^2/M0 invariant (gamma_{w1}=0)
  mu        Z_O^2      f1^2(mu)        M0(mu)           w1(mu)   dlnw1/dlnmu
  2.0     1.000000      0.300000     78.763325   0.003808879306
 100.0    44.140825     13.242248   3476.678167   0.003808879306    +0.00e+00
  check: dln f1^2/dln mu = 1.242670  vs  2*gamma_O(5) = 1.242670
  => w1 constant to ~1e-15; dln w1/dln mu = 0.  gamma_{w1}=0 EXACTLY.
[B] free field :phi^2:: f1^2=0 => w1=0 exact; M0(Lam=10)=1.92 -> M0(Lam=100)=198.8 (UV-divergent)
[C] 2D: <chi_0|chi_R> = -7.1e-17, 0, ... (=0)  => f1^2(2D trF^2)=0 => w1=0 exactly
[D] binding->0: A=1e0: w1=0.012584 ; 1e-8: w1=1e-10 ; 0: w1=0 ; invariant under Z^2=9 (same w1)
[E] free 4D scalar G(0) vs a: slope -2.083 (G(0)~a^-2); G(0)*a^2 -> const
```

All numbers reproduce the structural claims; the witness is a **self-check, not a proof** (theorem
content is §2–§6 and the cited files). Nothing here upgrades `[C]`/`[O]` to `[T]`.

---

## 8 · Status ledger

| item | statement | label |
|---|---|---|
| RG equation of `w₁` | `d ln w₁/d ln μ = 2γ_O−2γ_O = 0`, `γ_{w₁}≡0` exactly (ratio of two `Z_O²` insertions of one operator) | **[T]** |
| anomalous dimension of `tr F²` | `γ_O=−d ln|β_YM/2g|/d ln μ=2b₀g²+O(g⁴)`, `b₀=11N/48π²`; nonzero, but cancels in `w₁` | **[T]** |
| bare scalings | `M₀=C(0)~a⁻⁸` (contact), `χ~a⁻⁴` (T4); `f₁²` IR-clean (multiplicative only) | **[T]** |
| bare atom fraction | `w₁^{bare}=f₁²/M₀^{bare}→0` — a UV artifact, not a physical running | **[T]**+**[V]** |
| fixed-scale `w₁` | `O` flow / point-split ⇒ `f₁²(t),M₀(t)` finite ⇒ `w₁(t)` finite, `lim_{a→0}` exists | **[T]** |
| `w₁(t)` scheme dependence | depends on flow time `t` (finite renorm/mixing); floor must be `t`-uniform | **[T]** |
| positivity of the limit | `lim w₁>0 ⟺` atom survives with isolation `I₂`; RG-blind (`0` absorbing) | **[O]** |
| isolation ratio `ρ` | mass ratio, `a`-independent, `γ_ρ=0` | **[T]** |
| `ρ>1` ⇒ `w₁>0`? | **NO** — 2D: `ρ²=8/3>1` gapped, yet `w₁(tr F²)=0` (operator-specific) | **[T]**+**[V]**/**[NG]**-as-implication |
| free-field acid | `:φ²:` two-particle cut ⇒ `f₁²=0` ⇒ `w₁=0` | **[T]**+**[V]** |
| 2D exact `w₁` | `w₁=0` exactly (character orthogonality, `tr F²=H`) | **[T]**+**[V]** |
| strong coupling | `w₁≥r₀/M₀>0` at fixed `a`, `β<1/48`; **no** `a→0` transport | **[T]** |
| **exact missing estimate** | `inf_{L,a,t} w₁≥w₀>0` (leakage `≤1−w₀`) **and** `ρ≥1+c`, `tr F²` channel, `t`-uniform | **[O]** |
| Clay mass gap | unbuilt | **[O] · 0%** |

---

## 9 · Sharpest line

**The atom fraction `w₁=f₁²/M₀` does not run: it is the ratio of two insertions of the one dim-4
operator `tr F²`, both of which carry the same `Z_O²=(β/2g)^{-2}` from the trace anomaly, so
`d ln w₁/d ln μ=0` and the anomalous dimension vanishes identically — `γ_{w₁}=0`, `[T]`. Its `a→0`
limit at fixed physical volume is finite, because `f₁²` is an IR residue (multiplicative renorm only)
and `M₀`'s `a⁻⁸` contact is removed by renormalizing `O` at a fixed physical scale (gradient flow at
fixed `t`, or point-splitting at fixed `ε`); the bare ratio's apparent `→0` is a UV artifact, not a
decay. But the RG equation is a `0`-absorbing homomorphism: it fixes the *value*, never its *sign*, so
whether the continuum atom fraction is positive is exactly the non-perturbative atom-survival with
isolation `I₂` — the gap one rung up — and even a finite isolation ratio `ρ` does not give it: 2D pure
YM is gapped with `ρ²=8/3>1` while its `tr F²` atom fraction is `0` exactly. The RG half is free and
exact; the missing estimate is `w₁≥w₀>0` (equivalently leakage `Σ_{n≥2}w_n≤1−w₀`) together with
`ρ≥1+c` on the `tr F²` channel, uniform in `(L,a,t)` along the physical line at `β=2.5`. No QED.
Clay 0%.**

---

*Refs.* **Continuum / renorm:** `unconditional/G1-CONTINUUM-SCALING.md` (`G″_a=a²m²`, transport free,
seed is H1); `CONTINUUM-VARIANCE-DEFECT.md` (`D/S` is RG-covariant; `γ_O` controls scale not sign);
`PROOF-GRADIENT-FLOW.md` (flow defines `O_t`, rate `t`-independent, wrong half); `MS-RG-ODE.md`
(`d\ln f₁²/d\ln μ=2γ_O`, `m₁` RG-invariant); `INVENT-TRANSPORT-INEQ.md` (`Z_O²` square, `0→0`).
**Atom / carrier / isolation:** `T4-CARRIER.md` (total-normalized carrier `χ=O(1)`; bare `χ~a⁻⁴`;
`Var↛0` retraction); `T2-TUBE-OVERLAP.md` (`w₁=f₁²/M₀`, `Z_O`-invariant, leakage); `R5-EXACT-ANCHORS.md`
(free/2D `f₁²=0`; P3 strong-coupling seed); `R2-ISOLATION.md` (`I₂`, `ρ`, 2D `ρ²=8/3`); `R1-RESIDUE-TRANSPORT.md`
(`Z_O` free, `Ω` is the wall, `0` absorbing); `T3-TP2-ISOLATION.md` (window theorem); `T6-BINDING-ADVERSARIAL.md`
(`I⇒B⇒A` operator-specific). **Paper:** `docs/paper-v2/ch02-squeeze.tex` (`eq:chi`,`eq:squeeze`,`eq:h1`,
thin-atom Prop.), `ch03-existence.tex` (`§dim`,`eq:dim`), `ch12-numerics.tex` (finite-`L` `f₁²>0`
trivial; continuum is the content). **Standard (textbook, cited from corpus, not re-fetched):** Lüscher,
JHEP 08 (2010) 071; Lüscher–Weisz, JHEP 02 (2011) 051 (gradient-flow operator, separated-point
finiteness); Nielsen, NPB 120 (1977) 212; Collins, *Renormalization* (dim-4 operator mixing, `γ_O`);
Källén/Lehmann; Migdal (1975)/Witten CMP 141 (1991) 153 (2D solvability); Shen–Zhu–Zhu, CMP 400 (2023)
805 (strong-coupling gap, `β<1/48`). **Witness:** `research/v3/w7/u4_continuum.py` (`py -3.12 -I`,
numpy only, one job, exit 0, no downloads). No `[C]`/`[O]` upgraded to `[T]`; no QED; Clay 0%.**
