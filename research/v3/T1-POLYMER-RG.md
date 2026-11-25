# T1-POLYMER-RG — can a rigorous polymer / block-spin RG push the binding expansion to `β=2.5`? **No**, and the exact missing estimate is the **block effective activity** at the physical scale (equivalently the decay of connected `tr F²` polymer weights), missed by ~40–55×; the wall is the crossover, not a constant

> **2026-11-24 · Wave T1 (binding criterion, next species after R).** Attack: the atom `f₁²>0` is a
> **binding** phenomenon; the only exact anchor is P3 at `β<1/48`; the missing object is a binding
> theorem at `β=2.5`. **Question:** can the convergent strong-coupling (Kotecký–Preiss) polymer
> expansion, or a Gawedzki–Kupiainen-type block-spin/recentering reorganization, be **pushed to
> `β=2.5`**? Compute the minimal connected `tr F²` polymer weight and the Ursell ratio at `β=2.5`,
> estimate the effective block activity, and state the exact missing estimate.
> **Reads.** `v3/R-WAVE-SYNTHESIS.md`, `R5-EXACT-ANCHORS.md`, `R2-ISOLATION.md`;
> `unconditional/PROOF-P3-UNCONDITIONAL.md`, `P3-FULLREP-KP.md`, `P3-FOURWALL-DERIVATION.md`,
> `P3-DRESSED-MASS.md`; `v3/C1-SZZ-NSCALING.md`, `C2-BALABAN-UNIFORM.md`; `Q6-H1-STATUS.md`;
> `docs/paper-v2/ch05-theorems.tex` (`sec:p3`).
> **Witness.** `research/v3/w7/t1_polymer_rg.py` (`py -3.12 -I`, numpy+scipy.special.iv,
> deterministic, ~1 s, exit 0; output `t1_polymer_rg.out`; numbers quoted §1–§5). No third-party
> downloads.
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional · **[NG]** no-go for a
> named class · **[O]** exactly the missing step · **[V]** numerical witness. **No QED. Clay 0%.**

---

## §0 · Verdict (up front; no QED)

1. **No rigorous polymer expansion reaches `β=2.5` `[NG]` for the named class** (bare Kotecký–Preiss /
   Osterwalder–Seiler strong-coupling, and its single-layer block-spin/one-loop reorganization). The
   obstruction is not a bad constant — it is that **at `β=2.5` the strong-coupling expansion has no
   small parameter at all.** Three independent, mutually reinforcing diagnostics, all computed exactly
   for SU(2) [V]:

   | diagnostic | `β=1/48` (wall) | `β=1` | `β=2.5` | needs |
   |---|---|---|---|---|
   | trivial-character weight fraction `f₀=c₀/e^β` | **0.979** | 0.416 | **0.165** | `→1` (dilute gas) |
   | fundamental activity `u=c₁/c₀` | 0.0104 | 0.480 | **1.014** | `<1` (any KP) |
   | KP/Ursell ratio `λ_aΦ(1,β)` (d=4) | 3.14 | 419 | **4.6×10³** | `<1` (convergence) |

   At `β=2.5` only **17 %** of the plaquette weight is on the trivial irrep (the measure is *dense*
   with nontrivial flux), `u(2.5)≈1.01>1`, and the Ursell ratio exceeds its wall value by **1.5×10³**.
   The single-plaquette term alone already breaks the KP sum by two orders, **independent of the animal
   constant** (a calibration-free no-go, §1.4).

2. **The minimal connected `tr F²` polymer is still the four-wall tube, but at `β=2.5` its weight
   `u^{4T}` does not decay** [T]+[V]: `u(2.5)^{4}=1.06>1`, so the strong-coupling "mass"
   `m₁a=4\ln(1/u)` is `−0.057<0`. The expansion is not merely slowly convergent — it is in the **inverted**
   (dense, wrong-sign) regime. This kills any residue route through the tube series at `β=2.5`, not just
   the remainder (§2).

3. **A block-spin / one-loop reorganization does not create a small parameter** [NG]+[V]. Coarse-graining
   is asymptotically-free flow to the IR, `β_eff=β−8β₀\ln(ξ/a)` with `β₀=11N/(48π²)`; over the *physical*
   range `ξ/a≈1.3` this gives `β_eff=2.40`, `u_eff=0.988` — a `1.03×` drop. To reach the KP threshold
   `u<1/λ_a≈0.018` one would need `ξ/a≳760`; the lattice has `ξ/a≈1.3`. **Blocked by the absence of scale
   separation** (§3, §4): at `β=2.5`, `m_{0⁺⁺}a≈0.8`–`1.9`, so there are **0–1 RG steps** between the
   lattice and the correlation length.

4. **A Gawedzki–Kupiainen recentering is the *weak*-coupling face and fails from the other side** [NG].
   GK expands the fast fluctuation around the block vacuum; its small parameter is `1/β_eff≈0.4` at the
   physical scale — `O(1)`, not small. Additionally it computes the partition function, not the residue
   `f₁²`, and never reaches the crossover. So the two reorganizations are on **opposite sides** of
   `β≈2.5` and neither covers it (§4; the same "one species, two faces" as `C2 §4.2`).

5. **The exact missing estimate** [O] is a bound on the **effective block activity at the physical
   scale**, `u_eff(β=2.5)≤1/λ_a` — equivalently a **decay of connected `tr F²` polymer weights**,
   `|w(P)|≤e^{−κ|P|}` with `κ≥\ln λ_a` — i.e. `ρ_{Ur}^{eff}(2.5)<1`. Bare and one-loop values miss it by
   **~40–55×**, and closing it is **equivalent to controlling the crossover / the mass gap itself**; it is
   not a small-constant improvement and cannot be manufactured by any reorganization that keeps the
   expansion perturbative. This is the sharp, honest statement of the missing step (§6).

**One line.** *The bare strong-coupling expansion is dead at `β=2.5` (trivial fraction 17 %, activity
`≈1>1`, Ursell ratio `10³`); a block-spin/one-loop reorganization lowers the activity by only `1.03×`
because there is no scale separation (`ξ≈a`, 0–1 steps); GK recentering is the weak-coupling face and
equally fails; the exact missing estimate is the block effective activity at the physical scale,
`u_eff≤1/λ_a`, whose proof is the crossover/mass gap. No bare or single-layer block-spin/GK polymer
expansion reaches `β=2.5`. No QED. Clay 0%.*

---

## §1 · The expansion and its convergence criterion, computed exactly

### 1.1 Exact SU(2) character coefficients `[T]+[V]`

Wilson action `S=βΣ_p(1−½Re tr U_p)`, `e^{(β/2)Re tr U_p}=e^{β\cosθ_p}` (SU(2), `trU=2cosθ`). The
character expansion `e^{β\cosθ}=Σ_{n≥0}c_n(β)χ_{j}(θ)` with `n=2j=#boxes`, `dim=n+1`, has the exact
coefficients (Weyl/Itzykson–Zuber, `I_m` the modified Bessel function)

```
   c_n(β) = I_n(β) − I_{n+2}(β) ,     n = 0,1,2,…   (j = n/2) .
```

The completeness identity `Σ_n (n+1)c_n = e^β` is verified to `10^{-16}` at `β=1/48, 1, 2.5` (witness
[C]), so the coefficients are exact, not fitted. The **single-plaquette activity** is
`u(β):=c_1/c_0` (the fundamental, `n=1`), and the **per-plaquette KP activity at grading `a`** is

```
   Φ(a,β) = Σ_{n≥1} (c_n/c_0)(n+1) e^{a n} .
```

### 1.2 The Ursell / Kotecký–Preiss criterion `[T]`

Following `P3-FULLREP-KP §2`, the connected-polymer (Ursell) series for a correlation/free-energy
converges absolutely iff the Kotecký–Preiss sum through a fixed plaquette `x` is finite:
`S_x(a)=Σ_{P∋x}|w(P)|e^{a|P|}≤λ_aΦ(a)/(1−λ_aΦ(a))`, i.e. iff

```
   ρ_Ur(a,β) := λ_a Φ(a,β) < 1 ,        λ_a = O(Δ)  (animal/tree-growth constant, Δ = plaquette-adjacency degree) .
```

`Δ=4(2(d−1)−1)` (`=20,12,4` for `d=4,3,2`). The *literal* strong-coupling radius `1/48` is the
Shen–Zhu–Zhu / Bakry–Émery value (`C1 §1`, `β<1/(16(d−1))`); the crude KP constant `λ_a=eΔ` is tighter
(§1.4). Both are *method* radii; neither is a phase boundary. The point of §1.4 is that the no-go needs
**neither** — it is calibration-free.

### 1.3 The numbers `[V]`

Witness [E] (`λ_a=eΔ`, `a=1`):

```
   β        u=c₁/c₀   Φ(1)        ρ_Ur(d=4)   ρ_Ur(d=3)   ρ_Ur(d=2)
   0.0208   0.0104    5.8e−2      3.1         1.9         0.63      (the wall region)
   0.1      0.050     0.30        16          9.8         3.3
   1.0      0.480     7.7         4.2e2       2.5e2       84
   2.5      1.014     85          4.6e3       2.8e3       9.3e2
```

`ρ_Ur(2.5)/ρ_Ur(1/48)≈1.5×10³`. The threshold `u<1/λ_a=0.018` is reached only at
`β≈0.037`; `u(2.5)` is **55× over**. (The crude constant's own `ρ_Ur=1` crossing is at `β*≈0.0067`, so
`β=2.5` is `372×` beyond it in β. The SZZ wall `1/48` is looser; either way `β=2.5` is `10^1`–`10^3`
beyond.)

### 1.4 The calibration-free diagnostic — diluteness `[T]+[V]`

The KP criterion `λ_aΦ(1)<1` has `λ_a≥1` always, so a *necessary* condition is `Φ(1)<1`, hence (since
`e^{a n}≥1`) the total nontrivial activity `Φ(0)=Σ_{n≥1}(n+1)c_n/c_0 = (e^β−c_0)/c_0 < 1`. Equivalently the
**trivial-character weight fraction**

```
   f₀(β) := c_0(β)/e^β      must be close to 1  (a dilute gas of nontrivial plaquettes) .
```

Exact (witness [C]): `f₀ = 0.979 (β=1/48) → 0.416 (β=1) → 0.165 (β=2.5)`, and the nontrivial/trivial ratio
`Φ(0)=0.021 → 1.40 → 5.05`. **At `β=2.5` only 17 % of the weight is trivial** — the plaquette measure is
*dense* in colored flux. This is decisive **without any animal constant**: no per-plaquette
perturbative expansion can converge when the "small" (nontrivial) sector carries 83 % of the weight.
This is the cleanest single-number statement of the wall at `β=2.5`.

---

## §2 · The minimal connected `tr F²` polymer at `β=2.5` — and the sign inversion

### 2.1 The four-wall tube (restated, `P3-FOURWALL-DERIVATION §1`, `[T]`)

The minimal connected, gauge-invariant, nonvanishing `tr F²` polymer linking the `t=0` and `t=t`
insertions has `|P|≥4T` (`T=t/a`): each of the source plaquette's four boundary links must be contracted
(a shared partner plaquette), no partner covers two links, and only temporal plaquettes advance time, so
each of the 4 flux strands crosses `T` slices. The unique minimizers are the **straight four-walled
tubes**, weight

```
   w(B_T) = a₀(β)\,u^{4T} = a₀ e^{−m₁t},     m₁ a = 4\ln(1/u) ,   a₀>0  (positive recoupling) .
```

The next class is the perimeter-6 (1×2) tube, `w∝u^{6T}=u^{2}×w(B_T)` per slice; the two-tube/two-particle
class is `u^{8T}`. This is `P3`'s binding mechanism.

### 2.2 At `β=2.5` the tube weight is `u^{4T}` with `u>1` — inversion `[T]+[V]`

Witness [D]: `u(2.5)=1.014`, so `u^{4}=1.059>1` and `m₁a=4\ln(1/u)=−0.057<0`. **The minimal tube does not
decay.** The "binding" hierarchy of P3 — tube decays as `e^{−m₁t}`, competitors as `e^{−m₂t}` with
`m₂>m₁` — requires `0<u<1`; at `β=2.5` it is inverted, every additional plaquette *increases* the weight
until the animal entropy is cut off. So:

> **It is not only the *remainder* that is uncontrolled at `β=2.5`; the *leading* term has the wrong
> sign of decay.** A residue lower bound `f₁²=a₀(1−O(u))` is impossible to extract because `a₀u^{4T}` is
> not a decaying exponential, and `O(u)` is `O(1)`. `[T]+[V]`

(The small-`β` closing `m₁a=−4\ln(β/4)` of `P3-FOURWALL §4` extrapolates to `+1.88` at `β=2.5` only
because it *replaces* the exact `u(2.5)=1.01` by the small-β proxy `u≈β/4=0.625`; that is out of its
domain. The exact activity is the operative one.)

This is the sharpest structural statement of the wave: **binding is a small-β (dilute-tube) phenomenon;
at `β=2.5` the tube gas is dense and the perturbative binding picture does not merely fail to converge —
it signals a different, nonperturbative phase.**

---

## §3 · Block-spin / one-loop reorganization: does coarse-graining make the activity small?

### 3.1 The estimate `[T]+[V]`

Asymptotic freedom: `dg²/d\ln μ=−2β₀g⁴`, `β₀=11N/(48π²)`. With `β=4/g²` (SU(2)),
`dβ/d\ln μ=8β₀>0`. Coarsening a lattice to the physical IR scale `ξ` (flow `a^{-1}→ξ^{-1}`) gives

```
   β_eff(ξ) = β − 8β₀ \ln(ξ/a) ,      8β₀(SU(2)) = 0.372 .
```

Witness [F] (physical `ξ/a≈1.3`, §4):

```
   xi/a    β_eff    u_eff
   1.0     2.50     1.014
   1.3     2.40     0.988    ← physical range
   2.0     2.24     0.942
   4.0     1.99     0.861
   20      1.39     0.644
   800     0.017    0.008    ← where KP threshold is met
```

Over the physical range the activity drops `1.014→0.988` (**1.03×**). Reaching `u<1/λ_a=0.018` needs
`ξ/a≳760`; the lattice has `ξ/a≈1.3`. **Coarse-graining does not create a small parameter.**

### 3.2 Why (the flow is *tangent* to the crossover) `[T]`

Two ways to see it, both fatal:
- **AF flows the wrong way for "reach a controlled regime by coarsening."** Coarsening moves `β` *down*
  (toward strong coupling) — the eventual strong-coupling regime is at the *confinement* scale, which
  would need `ξ/a~10²–10³` to reach from `β=2.5`. At `β=2.5` the theory has `ξ/a~1.3`, so the flow never
  gets there; and flowing the *other* way (to the UV, larger `β`) is the weak-coupling / continuum
  direction, not the gap.
- **`β=2.5` is the crossover `β_c≈2.3`.** It sits where neither phase's expansion applies: not dilute
  strong coupling (`f₀=0.17`), not weak coupling (`1/β_eff≈0.4`). The RG map has no small parameter to
  deliver at the fixed point / crossover.

This is the constructive-RG face of `C2 §1.2`: the IR end of the tower is the gap itself (Channel I), and
`β=2.5` *is* the IR end.

---

## §4 · Gawedzki–Kupiainen recentering / two-scale — the weak-coupling face, and scale separation

### 4.1 GK recentering `[NG]`

GK write `U=\bar U·e^{iA}` (block vacuum × fast fluctuation) and expand the fast part `A`. The small
parameter is the fluctuation amplitude, controlled by the **stiffness** `β_eff`; convergence needs
`1/β_eff≪1`. At the physical scale `β_eff≈2.4`, `1/β_eff≈0.42` — `O(1)`. Moreover GK-type expansions
prove the **partition function / UV stability** (the `C2`/Balaban `L3` object), not the residue `f₁²`,
and they live at fixed volume at *weak* coupling — the opposite end from the confined IR where the atom
is bound. So GK recentering fails to reach `β=2.5` for the same reason the weak-coupling/perturbative
sector does: `β=2.5` is not weak enough.

### 4.2 No scale separation at `β=2.5` `[T]+[V]`

Witness [G]: the `0⁺⁺` glueball mass in lattice units at `β=2.5`, two estimates —
strong-coupling leading `m a = 4\ln(4/β)=1.88` (overestimates at the crossover) and scaling-region
`m_{0⁺⁺}/\sqrtσ≈3.5`, `\sqrtσ a≈0.22 ⇒ ma≈0.77` — give `ξ/a = 1/m a ∈[0.53,1.30]`. So

```
   number of block-spin (b=2) steps between a and ξ  =  \log_2(ξ/a) ∈ [−0.9, 0.4]  ⇒  0–1 steps .
```

**A block-spin RG has essentially no tower to run at `β=2.5`.** The "two-scale" expansion's two scales
coincide. This is the geometric reason the reorganization cannot help: there is no intermediate scale on
which to define a controlled effective theory.

---

## §5 · Tests on solvable / limiting cases

### 5.1 2D pure Yang–Mills (solvable) — reaches all `β`, atom `=0` `[T]`

2D pure YM is exactly solvable (Migdal/Witten) and **confining for every `β`** with an exact area law —
it "reaches `β=2.5`" trivially, by solvability, not by any polymer expansion. Yet its `tr F²` atom is
**`f₁²=0` exactly** (`R5 §1.1`: `tr F²=tr E²=H`, conserved and diagonal; `⟨0|trF²|G⟩=C₂(G)⟨χ_0|χ_G⟩=0`).
Two lessons, both favorable to this file:
- The diagnostic agrees the *strong-coupling* expansion is dead at `β=2.5` in `d=2` too
  (`ρ_Ur=9.3×10²`, witness [E]); 2D is solved by the *Migdal character* route, not the KP
  cluster expansion.
- **Reaching `β=2.5` does not give the atom**: 2D reaches it and `f₁²=0`. So "pushing an expansion to
  `β=2.5`" is at best *necessary*, never sufficient, for binding — exactly `R5`'s ceiling/atom split.

### 5.2 3D pure Yang–Mills (super-renormalizable) — same failure `[T-standard]+[V]`

The activity is a per-plaquette SU(2) object independent of `d`; only `Δ` changes (`Δ=12`). Witness [E]:
`ρ_Ur(2.5,d=3)=2.8×10³`, `u(2.5)` identical `1.01`. So the bare expansion fails at `β=2.5` in 3D too.
Since 3D is *super-renormalizable* (a strictly "easier" theory: trivial continuum, dimensionful coupling,
continuum limit `β→∞`), this shows the `β=2.5` failure is **dimensional-robust**, not a 4D accident — a
moderate-coupling lattice at `β=2.5` is beyond any bare expansion in every `d≥2`.

### 5.3 Free-field acid test — must fail; it "fails correctly" `[T]+[V]`

Witness [H]: as `β→∞`, `u=c₁/c₀→2` (both `I_n~e^β/\sqrt{2πβ}`) and `Φ(1)→∞`, so `ρ_Ur→∞`
(`β=10`: `ρ_Ur=5.8×10⁶`; `β=50`: `7.8×10¹⁷`). The activity diagnostic correctly reports **"no expansion"**
on the free field — the strong-coupling expansion is maximally non-convergent there. And the free
massless composite `tr F²` has `f₁²=0` structurally (two-gluon cut, `R5 §1.6`). So the no-go **never
manufactures a false atom**: it returns "no expansion" where the true answer is "no atom," and the acid
test passes trivially (the tool is a *no-go*, not an atom-claimer).

### 5.4 What the tests establish

| test | true answer | diagnostic | verdict |
|---|---|---|---|
| 2D YM `β=2.5` | solvable, `f₁²=0` | expansion dead (`ρ_Ur~10³`) | pass — no false atom; reach≠atom |
| 3D YM `β=2.5` | confining, gap | expansion dead (`ρ_Ur~10³`) | pass — dimensional-robust failure |
| free field `β→∞` | `f₁²=0` | expansion maximally dead (`ρ_Ur→∞`) | pass — acid (no false atom) |

---

## §6 · The exact missing estimate (the deliverable)

### 6.1 Statement `[O]`

> **`(M1)` — effective-activity / polymer-weight decay at the physical scale.** There exist
> `λ_a>0` (the animal constant) and `β`-independent constants such that, **uniformly through the
> physical window including `β=2.5`**, either
> ```
>   (polymer-weight form)   |w(P)| ≤ e^{−κ|P|}   for every connected tr F² polymer P,  with  κ ≥ \ln λ_a ;
>   (block-activity form)   u_eff(β=2.5) ≤ 1/λ_a        (the coarse-scale activity at scale ξ) ;
>   (Ursell form)           ρ_Ur^{eff}(β=2.5) := λ_a Φ(u_eff) < 1 .
> ```

All three are the same statement. Bare values: `κ_{bare}(2.5)=−\ln u≈−0.014<0` (the weight *grows*);
`u(2.5)=1.01` vs the needed `1/λ_a≈0.018` (**55× over**); one-loop improvement to `u_eff=0.988`
(**still 54× over**). So `(M1)` is missed by `≈40–55×` in the activity, and by an *unbounded* factor in
the Ursell ratio (`4.6×10³` vs `<1`).

### 6.2 Why it is the *same* missing estimate as the gap (the honesty clause) `[NG]`

`(M1)` is **not** a small-constant improvement any reorganization can supply. Three reasons:
1. **It is calibration-free necessary** (`§1.4`): `f₀(2.5)=0.165≪1` — the weight is dense in colored
   flux; no per-plaquette perturbative expansion converges there, independent of `λ_a`.
2. **The leading term is inverted** (`§2.2`): `u(2.5)>1` — the tube does not decay, so `(M1)`'s `κ>0`
   cannot come from the dilute-tube picture at all; it must come from a genuinely nonperturbative
   mechanism that *makes* the effective activity collapse.
3. **That mechanism is confinement/binding** — i.e. the mass gap. A block-spin flow that drives
   `u_eff→0` *is* the statement that the theory is in the confining phase with a mass gap; proving it is
   the crossover control (`C2 §4`: the tower's IR end = H1 = the gap). So `(M1)` at `β=2.5` is
   **equivalent to the gap, not weaker** — supplying it perturbatively is circular.

Equivalently, in the language of the R wave: the "binding theorem at `β=2.5`" that `R-WAVE-SYNTHESIS §7`
names as the next species **cannot be a Kotecký–Preiss/block-spin statement at all** — those are
small-coupling tools on both sides of the crossover, and `β=2.5` is the crossover. A genuine binding
theorem must be a *non-perturbative* statement (variational/spectral/order-parameter), not an expansion;
the expansion is exactly what `(M1)` says does not exist.

### 6.3 What would and would not close it (constructive, but honest) `[O]`

- **Would NOT:** optimizing `λ_a`, choosing a better animal bound, a `b=2` block-spin step, or a GK
  recentering at `β_eff≈2.4`. Each moves `u_eff` by `O(10\%)` into a `40–55×` gap (`§3–§4`). The failure
  is structural (dense gas + no scale separation), not quantitative.
- **Would (a target, not a result):** a **non-perturbative effective-activity contraction** — a proof
  that the physical-scale measure has `f₀→1` (dilute) as a *consequence of confinement* at `β=2.5`, e.g.
  from an area-law / center-vortex / monopole-condensate order parameter that pins the colored sector.
  This is a *different tier* (order parameter, not expansion), and it is exactly what the R wave said the
  variance/measure tools cannot reach. It is `[O]`, unbuilt; recording it as the precise target, not a
  result.

---

## §7 · Status ledger

| item | statement | label |
|---|---|---|
| exact SU(2) coefficients | `c_n=I_n−I_{n+2}`, `Σ(n+1)c_n=e^β` (verified `10^{-16}`) | **[T]+[V]** |
| four-wall tube weight | `w(B_T)=a₀u^{4T}`, `m₁a=4\ln(1/u)` (P3 geometry) | **[T]** |
| activity at `β=2.5` | `u(2.5)=1.014`, `u^4=1.06>1`, `m₁a=−0.057<0` (train inverts) | **[T]+[V]** |
| diluteness at `β=2.5` | `f₀=c₀/e^β=0.165` (17 % trivial; dense gas) | **[T]+[V]** |
| Ursell/KP ratio | `ρ_Ur(2.5)=4.6×10³` (d=4), `1.5×10³×` the wall value | **[V]** |
| bare KP reaches `β=2.5`? | **No** — necessary condition `f₀→1`, `u<1` fails | **[NG]** |
| block-spin / one-loop RG | `β_eff(2.5;ξ/a=1.3)=2.40`, `u_eff=0.988`; needs `ξ/a≳760` | **[NG]+[V]** |
| GK recentering | weak-coupling face, `1/β_eff≈0.4` not small; gives `Z`, not `f₁²` | **[NG]** |
| scale separation `β=2.5` | `ξ/a∈[0.53,1.30]`, 0–1 RG steps — no tower | **[T]+[V]** |
| 2D test | solvable all `β`, `f₁²=0`; expansion dead at 2.5 — reach ≠ atom | **[T]+[V]** |
| 3D test | `ρ_Ur(2.5,d=3)=2.8×10³` — dimensional-robust failure | **[T-standard]+[V]** |
| free-field acid | `u→2`, `ρ_Ur→∞`; `f₁²(free)=0` — no false atom | **[T]+[V]** |
| **`(M1)` missing estimate** | `u_eff(2.5)≤1/λ_a` / `κ≥\ln λ_a` / `ρ_Ur^{eff}<1` — missed `40–55×` | **[O]** |
| `(M1)` vs the gap | equivalent, not weaker; perturbative supply is circular | **[NG]** |
| polymer/block-spin RG to `β=2.5` | **does not exist** | **[NG]** |
| Clay mass gap | unbuilt | **[O] · 0%** |

---

## §8 · Sharpest line

**No rigorous bare or single-layer block-spin/GK polymer expansion reaches `β=2.5`, and the reason is the crossover, not a
constant.** The exact SU(2) activity `u(β)=c₁/c₀` is `0.0104` at the P3 wall and `1.01` at `β=2.5`; the
trivial-character weight fraction falls `0.98→0.17`; the Ursell ratio `λ_aΦ(1,β)` rises `3.1→4.6×10³`.
The minimal `tr F²` polymer is still the four-wall tube, but at `β=2.5` its weight `u^{4T}` does not decay
(`u⁴=1.06`, `m₁a<0`) — the dilute-tube binding picture inverts. A block-spin/one-loop reorganization
lowers the activity only `1.03×` over the physical range because there is **no scale separation**
(`ξ/a≈1.3`, 0–1 RG steps), and GK recentering is the *weak*-coupling face (`1/β_eff≈0.4`), equally
outside `β=2.5`. Tests are clean: 2D (solvable all `β`, `f₁²=0`) and 3D (super-renormalizable) both show
the expansion dead at `β=2.5`, and the free field (`u→2`, `f₁²=0`) passes the acid test with no false
atom. The **exact missing estimate** is the physical-scale effective activity, `u_eff(2.5)≤1/λ_a`
(equivalently `|w(P)|≤e^{−κ|P|}`, `κ≥\ln λ_a`; `ρ_Ur^{eff}<1`) — missed by `≈40–55×`, and **equivalent
to the mass gap**, so no perturbative reorganization can supply it. **The binding theorem the R wave asks
for cannot be an expansion at `β=2.5`; it must be non-perturbative. No QED. Clay 0%.**

---

## §9 · Witness — `research/v3/w7/t1_polymer_rg.py`

`py -3.12 -I t1_polymer_rg.py` (exit 0; numpy 2.4.6 + scipy 1.18.1 `special.iv`; deterministic, no MC,
~1 s; saves `t1_polymer_rg.out`). Quoted output (abridged):

```
[C] beta=2.5: SUM(n+1)c_n=12.18249396=e^beta (rel.err 9e-16); f0=c0/e^beta: 0.979(1/48) 0.416(1) 0.165(2.5)
[D] beta=2.5: u=1.01439, u^4=1.0588, m1a=-0.05715  ;  beta=1/48: u=0.01042, u^4=1.2e-8, m1a=18.26
[E] rho_Ur(1,beta)=lambda_a*Phi: beta=1/48: 3.14/1.89/0.63 (d=4/3/2) ; beta=2.5: 4.63e3/2.78e3/9.26e2
    crude KP needs u<0.01839 ; u(2.5)=1.01439 => 55.1x over ; rho(2.5)/rho(1/48)=1473x
[F] beta_eff(2.5; xi/a=1.30)=2.4025, u_eff=0.98773 (1.03x drop) ; need xi/a>~758 ; physical xi/a~1.3
[G] ma=1.880 (strong-cpl) / 0.770 (scaling) => xi/a=0.532/1.299 ; log2(xi/a)~0.38 => 0-1 RG steps
[H] free field: beta=10 u=1.708 rho_Ur=5.8e6 ; beta=50 u=1.940 rho_Ur=7.8e17 ; u->2, f1^2(free)=0
```

Every number reproduces the analysis of §1–§5. The witness is a **self-check, not a proof**: all theorem
content (the four-wall minimality, the KP criterion, the AF flow, the analyticity wall) is in §1–§6 and
its cited P3/C2 files. `[V]`

*Refs.* `unconditional/PROOF-P3-UNCONDITIONAL.md` (residue `f₁²≥r₀(β)`, `β<1/48`, `[T]` mod OS/KL);
`P3-FULLREP-KP.md` (`λ_aΦ<1`, `K≤CN`, full-rep KP); `P3-FOURWALL-DERIVATION.md` (four-wall tube,
`m₁=4\ln(1/u)`); `P3-DRESSED-MASS.md` (decoration resummation); `v3/C1-SZZ-NSCALING.md` (SZZ radius
`1/48`, `K_S∝N`); `v3/C2-BALABAN-UNIFORM.md` (Channel I = IR tower bottom = the gap; `(★)`; Lemma B);
`v3/R-WAVE-SYNTHESIS.md §7` (next species = binding criterion, not another variance bound);
`v3/R5-EXACT-ANCHORS.md` (binding is the discriminant; 2D/`β=0`/free are atomless);
`v3/R2-ISOLATION.md` (`I2` ≡ gap one rung up). **Standard inputs (textbook, not re-fetched):**
Kotecký–Preiss, CMP **103** (1986) 491; Osterwalder–Seiler, Ann. Phys. **110** (1978) 440;
Seiler, LNP **159** (1982) §II; Münster, PLB **95** (1980) 59, NPB **190** (1981) 439;
Gawedzki–Kupiainen, CMP **99** (1985) 197 / JSTAT (1984); Gross–Mezard / Itzykson–Zuber (SU(2)
character coefficients); Migdal (1975), Witten CMP **141** (1991) 153 (2D solvability).
**No third-party sources downloaded; witness run once, deterministic; no `[C]` upgraded to `[T]`;
no QED; Clay 0%.**
