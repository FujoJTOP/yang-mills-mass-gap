# INVENT-PERCOLATION — the random-cluster functional IS a genuine new species (it carries FK–FKG, which Peierls lacks), but its threshold is the *confinement* (1-form) transition and it needs a `q`-state that full-rep `SU(N)` does not have — so it dies at the same three walls

> **VERDICT (no QED).** A Fortuin–Kasteleyn / random-cluster percolation functional built on the gauge
> measure is a **genuine new species**, not a Peierls relabelling: it satisfies the exact identity
> **`Q_p(0↔x) = ⟨σ_0σ_x⟩`** (connectivity = correlator, verified to `10⁻¹⁴`) and it restores
> **FK–FKG positive association** — the one monotonicity the gauge measure lacks. **But** (i) its
> percolation threshold is the **order–disorder transition = the *confinement* (1-form) transition**, and at
> zero `T` 4D pure YM that threshold is **absent** (single phase, `B2`) ⇒ no threshold to measure the `0⁺⁺`
> gap from; (ii) the `q`-state FK representation **does not exist for full-rep `SU(N)`** (character weight ≠
> a `q`-state) ⇒ no FK–FKG there; (iii) where it exists (`Z_N`), the connectivity decay rate **is** the
> surface tension = the Peierls `σ_w` threshold. **A real new tool that dies at the same three walls.**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-num]** new numeric · **[V-here]** new this round.
Iron laws: **IL-i** acid (free massless `=0`, `:φ²:` discriminatingly `=0`); **IL-ii** nonabelian; **IL-iii** no small parameter at `α_s~O(1)`.

---

## §1 · What was built (the FK / random-cluster functional) `[T]`

On a lattice with edge set `E`, the Fortuin–Kasteleyn random-cluster measure (cluster weight `q`,
edge weight `p`) is
```
Z_RC(p,q) = Σ_{A⊆E} p^{|A|} (1−p)^{|E|−|A|} q^{k(A)},     k(A) = #connected components.
```
For `q=2` it is the Ising model (`p = 1−e^{−2β}`), and the **central object** is the connectivity
```
Q_p(0↔x) := Prob_p(0 and x lie in the same FK cluster).
```
**Exact identity (witness, 3×3 lattice, brute force, diff `≤1.2×10⁻¹⁴`):**
```
Q_p(0↔x) = ⟨σ_0 σ_x⟩.
```
So the "`0⁺⁺` cluster connectivity" **is** a correlator — and its **decay rate is the spectral gap**. Three
facts verified `[V-num]`:
- **(A)** `Z_RC(p,2) = Z_Ising(β)` up to a β-dependent prefactor (standard FK relation);
- **(B)** `Q_p(0↔x) = ⟨σ_0σ_x⟩` exactly (the point of the whole construction);
- **(C)** `Q_p(0↔x)` is **non-decreasing in `p`** (FK–FKG / positive association) — verified monotone.

**This is genuinely NOT a Peierls relabelling.** Peierls is a contour/surface-tension sum; percolation is a
*connectivity / correlation-length* functional, and — crucially — FK–FKG supplies a **monotone coupling**
that the contour argument never had.

---

## §2 · Why it does not lay the roof: three exact walls

### (i) `[NG]` channel: the threshold is the *confinement* (1-form) transition, and at zero `T` it is absent.
The connectivity `Q_p(0↔x)` is the **magnetization** correlator. Its percolation transition (`p_c`, where
`ξ_p→∞`) is the **order–disorder transition**. For the gauge theory the analogous percolating object is the
**center/magnetic (1-form)** sector, and its transition is **confinement/deconfinement** — *not* the `0⁺⁺`
gap. At
zero `T`, 4D pure YM is a **single phase** (`ROOF-TRANSFER` B2, `[T]`): there is **no confinement transition
to cross**, hence **no percolation threshold**, hence **nothing to measure the `0⁺⁺` gap as a distance
from**. The witness shows exactly this in the Ising toy: `Q→1` (percolates) for `β>β_c` — the
order-parameter behavior — while the `0⁺⁺`/energy gap is finite on *both* sides and vanishes only *at* the
(absent, in 4D zero-`T`) critical point.

### (ii) `[NG]` no `q`-state: full-rep `SU(N)` has no random-cluster representation.
`Z_RC` uses a cluster weight `q^{k(A)}` with **integer (or half-integer) `q`** — a `q`-state Potts/Ising
structure. The full-rep `SU(N)` Wilson weight is a **character sum** `Σ_r d_r χ_r(U)`, not a `q`-state
spin weight; there is **no FK expansion**, hence **no FK–FKG**, for full-rep `SU(N)`. The one genuinely new
ingredient (cluster-level positive association) is **exactly what is unavailable** — the same wall that
kills FKG in `BASE-FIELDMEASURE` R2.

### (iii) `[O]` where it exists (`Z_N`), the rate IS the Peierls threshold.
For `Z_N` gauge the FK representation exists, but the connectivity decay rate `1/ξ_p` **equals** the
surface tension `σ_w` of the dual (vortex-surface) cluster — i.e. the same `c₀` threshold as
`ROOF-PEIERLS-0FORM`. So percolation gives the **same bound** `σ_w > c₀` (blocked at `σ(2.5)=0.57 <
c₀≈1.1–1.95`), not a new one.

---

## §3 · The one construction that IS a genuine `P1` species (and its exact gap) `[V-here]/[O]`

The `q`-state FK is the wrong object for full-rep `SU(N)`, but percolation suggests a **level-set / field-intensity
percolation** functional directly on the `0⁺⁺` operator:
```
define a "0⁺⁺ cluster" as a connected set of plaquettes with  tr F² > θ.
P(θ) := percolation probability of these clusters;
𝔉_perc(μ; θ) := the decay rate of the cluster-density correlator.
```
This is a **genuine field-measure functional** (`PROFILE` P1): nonlinear in `μ`, defined at **every** β
(no transition needed), phase-sensitive through `θ`, and **acid-testable** (on the free massless field
`tr F²` is a Gaussian fluctuation with no percolating clusters at any fixed `θ` ⇒ `𝔉_perc=0`; checking the
`:φ²:` discrimination is the real acid test). It is **not** a Markov-semigroup functional, so Theorem BE
does not kill it, and **not** an order parameter, so B2 does not kill it.

**But its exact gap `[O]`:** there is **no proven inequality** from `𝔉_perc(μ;θ)` to any of the three
names of the strut (`χ_lattice ≲ a²`, `σ_w(2.5)>c₀`, `R_vac` monotone). The level-set percolation
probability is a **new observable**, computable in principle, but the link "percolation ⇔ floored
susceptibility" is open — and the natural reason it is open is the same one everywhere: at `α_s~O(1)` there
is no small parameter controlling the level-set threshold `θ`. **Exact failing step: a percolation-to-`χ`
inequality** — i.e. the strut again, in a fourth coordinate.

---

## §4 · Verdict

- **Genuinely different from Peierls?** **Yes as a representation** — it carries the exact
  `connectivity = correlator` identity and FK–FKG positive association, which the contour argument lacks.
  **No as a bound** — the `q`-state version detects the 1-form/confinement transition (absent at zero `T`),
  fails to exist for full-rep `SU(N)`, and where it exists reproduces the `σ_w` threshold.
- **New species?** The **level-set percolation** functional `𝔉_perc` is a genuine `P1` species (non-Markov,
  non-order-parameter, phase-sensitive, acid-testable) — the first construction this round that clears
  those three paradigm no-gos. **But it supplies no bound**: the percolation-to-`χ` inequality is exactly
  the missing strut, in a fourth coordinate.
- **Exact failing step pinned:** (1) `q`-state FK absent for full-rep `SU(N)`; (2) percolation threshold =
  1-form/confinement, absent at zero `T` (B2); (3) no percolation-to-`χ` inequality for the level-set
  functional.

**Not a proof; not a dead end.** It is the first invented species that is **not** an order parameter, **not**
a Markov functional, and **not** topological — and it localizes the wall to a single missing inequality:
*percolation probability ⇒ floored susceptibility*. **Clay 0%; no QED.**

*Witness:* `w7/invent_percolation.py` (3×3 lattice, brute force; FK identity up to prefactor; connectivity
`=⟨σσ⟩` to `1e−14`; FK–FKG monotone; percolation in the ordered phase). *Refs:* `ROOF-PEIERLS-0FORM.md`,
`ROOF-TRANSFER.md` (B2), `BASE-FIELDMEASURE.md` (Theorem BE, FKG wall), `PROFILE.md` (P1).
