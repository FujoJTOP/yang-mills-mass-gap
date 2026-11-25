# Y4-FESHBACH — TEST of the **gauge-invariant nonlinear Feshbach / multiscale RG projection**: it is a clean, exactly gauge-independent, non-vacuous **reformulation**, but it does **not** improve on the WG estimator — its `k=0` effective action is **Gaussian in the thermodynamic limit** (CLT) and equals the WG variance; the only genuine increment (all-`k`) is the **standard** gauge-invariant effective mass (estimate, not a bound); on the Savvidy cell it **inherits** the physical tachyon. Kill criterion **(a) repackages WG/H1** and **(c) cannot control the collective mode** both **TRIGGERED**; **(b) needs a gauge choice** is **NOT** triggered.

> **2026-11-24 · Wave Y4 — evaluating a proposed *new* machine.**
> Proposal under test (task): construct a gauge-**independent** Feshbach / multiscale RG projection that
> tracks the low-frequency `k=0` collective mode **without gauge fixing**; does it give a non-vacuous
> carrier with `Var ≠ V^{-1}`, and does it **improve** on the WG estimator? Tests on 1D/2D solvable
> models and the Savvidy cell; acid test free field `→0`/known and 2D `trF²`; kill criteria: repackages
> WG/H1, needs a gauge choice, cannot control the collective mode.
> Reads: `H1-RG-DESCENT.md` (`G″_j=(P A_j^{-1}P^T)^{-1}`, `G″=E[Σ]−D`, `η_j` breakdown at the Savvidy
> soft mode), `S3-UV-FESHBACH.md` (`G″=E[Σ]−D`, sign of `D`), `T4-CARRIER.md` (carrier `Var=χ_O` vs
> `χ_O/V`; thin-atom no-go T4-4), `WG-INVARIANT-ESTIMATOR.md` / `WG-L12-L16.md` / `WG-FSS-FIT.md`
> (`χ_O=0.328,0.377,0.411,0.469,0.587` L=4,…,16; slope `0.37±0.07`), `docs/paper-v2/ch02-squeeze.tex`
> (`eq:chi`, `eq:squeeze`, `eq:ward` = `G″·Var=1` Gaussian-only, §`G″ vs 1/χ`, thin-atom Prop.).
> **Witness:** `research/v3/w7/y4_feshbach.py` (`py -3.12 -I`, numpy only, one process; output §7).
>
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional · **[NG]** no-go for a
> named class · **[O]** exact missing step · **[V]** numerical witness here. **No QED. Clay 0%.**

---

## 0 · Verdict (up front; no QED, Clay 0%)

The proposal is **well-defined and passes the two "hygiene" tests** — it is exactly gauge-independent
and non-vacuous — **but it fails the task's substantive test**: it does **not** improve on WG in any
load-bearing sense, and it does **not** control the Savvidy collective mode. Precisely:

1. **Gauge-independent, no gauge fixing `[T]`.** Take a gauge-invariant density `O` (`tr F²`,
   plaquette, APE-smeared) and its total-normalized zero mode `q_O = V^{-1/2}Σ_x(O−⟨O⟩)`. The Feshbach
   effective action `S_eff(q) = −log E_μ[δ(q−q_O)]` is a functional of a gauge-invariant random
   variable, hence exactly invariant under every gauge transform, **with no gauge choice** (§1, §7[F]
   `|ΔQ|/|Q| = 0` per config; the link-space log mode drifts `O(1)`). **Kill criterion (b) NOT triggered.**

2. **Non-vacuous carrier, `Var ≠ V^{-1}` `[T]` — but it is exactly the T4/WG carrier.** `Var(q_O) = χ_O
   = Σ_z C(z)`, an `O(1)`, `L`-independent quantity, not `χ_O/V`. This is literally the WG estimator's
   object (`T4-CARRIER`, `WG-INVARIANT-ESTIMATOR`). So the "carrier" is **not new**: the Feshbach adds
   no carrier beyond H1's own `q` (§2).

3. **Does it improve on WG? NO `[NG]`.** Three independent reasons:
   - **Gaussian limit = reciprocal WG `[T]`.** In the Gaussian case the Feshbach curvature is
     `G″=S_eff''(q*)=(P A^{-1}P^T)^{-1}` and `Var=P A^{-1}P^T`, so `G″·Var=I` **exactly** (§7[A], `1e-15`;
     ch02 `eq:ward`). The Feshbach curvature *is* the reciprocal susceptibility in the only case it is
     exact.
   - **Off-Gaussian split collapses by CLT `[T]/[V]`.** The curvature `G″=S_eff''` and the reciprocal
     variance `1/χ` differ at finite `V` (one-plaquette: `G″=1` vs `1/χ=4`, factor 4; free composite
     `:φ²:`: ratio `1.37→1.11→1.04` at L=16,32,64), **but the difference is a non-Gaussianity that
     vanishes**: the excess kurtosis of `q_O` is `κ₄=3.92→1.79→0.84→0` (CLT), and `S_eff''→1/χ`. The
     whole **nonlinear** Feshbach effective action is **Gaussian in the thermodynamic limit**,
     `S_eff(q)=q²/(2χ_O)+o(1)`; all cumulants `κ_n=O(V^{1−n/2})→0`. So the nonlinear machinery carries
     **no thermodynamic content beyond `χ_O`** — which is H1/WG (§3).
   - **Ceiling, not pole `[T]`.** For a composite (continuum KL measure) both `G″` and `1/χ` are
     ceiling-type and `≠ m₁²` (free `:φ²:`: `G″≈1/χ≈0.24–0.33` vs `(2m)²=1`); they equal the pole only
     for a single-atom (elementary) KL measure. Same as ch02's `G″ vs 1/χ` (measured 2–5% at β=2.5) and
     `χ ≠ pole`. **Kill criterion (a) TRIGGERED.**

4. **Does it control the collective mode? NO `[NG]`.** The gauge-invariant collective coordinate is
   `q ~ B²` (it tracks the background amplitude). On the **Savvidy one-loop cell** its effective
   potential is `V_eff(B²)` with `V_eff''(B)<0` for small `B` — the **inherited tachyon** (§7[E]: sign
   flip at `B_c`, Savvidy minimum at `B_min>B_c`). The projection does **not** flip/stabilise it; the
   full-measure positivity (`1/χ_O = 1.7–3.0 > 0`, WG) is a **measure** fact (ensemble average / compact
   Haar), not a projection fact. In particular the gauge-invariant Feshbach **sidesteps** the
   `η_j≥1` breakdown of `H1-RG-DESCENT` only by **giving up mode-by-mode curvature control** and
   returning the ceiling `1/χ_O`. **Kill criterion (c) TRIGGERED.**

5. **Multiscale = the standard effective mass `[T-standard]/[C]/[V]`.** Doing the projection at every
   `k` gives `C(k)=Var(q(k))=G̃(k)`, the full gauge-invariant propagator; the effective mass from the
   `C(x)` decay recovers `m` (§7[G]). This is the **standard** gauge-invariant effective mass
   (`INVENT-GAUGE-INVARIANT-MASS.md`), an **estimate** `[C]/[V]`, **not** a lower bound, and **not new**.
   It genuinely adds the `k≠0` slice that WG's single number `χ_O=C(0)` lacks — but that increment is the
   established propagator/effective-mass tool, and the **edge** (branch point of the analytic `C(k)`) is
   the R3 circular/infinite-data target, unchanged.

**One line.** *The gauge-invariant nonlinear Feshbach is a correct, exactly gauge-independent,
non-vacuous **reformulation**: at `k=0` its curvature is the reciprocal WG susceptibility in the
Gaussian limit and collapses to it by CLT off-Gaussian (the "nonlinear" part is `O(V^{-1/2})` noise with
no thermodynamic content); at all `k` it is the standard gauge-invariant effective mass (an estimate);
on the Savvidy cell it inherits the physical tachyon. It reproduces H1/WG as a repackaging and cannot
control the collective mode — the two kill criteria that matter. The atom `f₁²≥r₀` remains `[O]`.
No QED. Clay 0%.*

---

## 1 · The proposal, stated once `[T]`

**Object.** Fix a gauge-invariant `0⁺⁺` density `O` (`tr F²`, the plaquette energy density, an
APE/gradient-flow-smeared `tr F²` at a fixed physical scale — the `T4`/`WG` carrier). Its
**total-normalized zero mode** and the **Feshbach effective action** are
```
   q_O[U] = V^{-1/2} Σ_x ( O(x) − ⟨O⟩ ) ,        e^{−S_eff(q)} = E_μ[ δ( q − q_O[U] ) ] ,          (1.1)
```
i.e. `S_eff = −log` of the marginal density of `q_O`. This is exactly the **Feshbach / RG map
integrated to the IR** (all high-frequency modes integrated out, leaving the `k=0` collective
coordinate). The **multiscale** version does the same for every momentum mode,
`q_O(k) = V^{-1/2}Σ_x e^{−ikx}(O−⟨O⟩)`, giving the family `S_eff^{(k)}` with curvature
`C(k) = Var(q_O(k))`.

**Why it is gauge-independent (no gauge choice) `[T]`.** `O` is a pointwise gauge-invariant composite,
`V^{-1/2}` is a constant, and `⟨O⟩` is a `c`-number; hence `q_O` — and every functional of it, in
particular `S_eff(q)` and its derivatives — is invariant under every gauge transformation. There is
**no** projection onto a gauge orbit and **no** gauge-fixing. (Contrast the **link-space** Feshbach of
`H1-RG-DESCENT`/`S3`, whose `Hess V` lives on `log U` and is gauge-dependent; §7[F] shows the
observable `Q` invariant to `0` and the link-log mode drifting `O(1)`.)

**Two readings of "nonlinear".**
- *(i) Curvature.* `G″ := S_eff''(q*)` = the Hessian of the effective action at its mode — the object
  `H1-RG-DESCENT` calls `G″_j`. Below we test whether this improves on `1/χ`.
- *(ii) Full shape.* The entire function `S_eff(q)` (= all cumulants of `q_O`) — the genuinely
  "nonlinear" content. Below we test whether the shape carries edge information the variance lacks.

Any **nonlinear collective coordinate** `Q = f(volume-averaged gauge-invariant observables)` is still
a `k=0` thermodynamic object and is covered by the same argument (§3.4).

---

## 2 · Q1 — non-vacuous carrier `[T]`, but it is the T4/WG carrier (nothing new)

`q_O` is the total-normalized zero mode. By `T4-CARRIER` Thm T4-1 (and §1 of `WG-INVARIANT-ESTIMATOR`),
```
   Var(q_O) = (1/V) Σ_{x,y} C(x−y) = Σ_z C(z) = χ_O   (O(1), L-independent) ,
   Var(q̄_O) = χ_O / V                                 (mean-normalized, self-averaging → 0) .          (2.1)
```
So the Feshbach carrier **is non-vacuous, `Var=χ_O ≠ V^{-1}`** — the acid test `[V]` confirms it
(free massive scalar `Var=1/m²` exactly at every `L`; §7[B]). **But this is precisely the `T4`/`WG`
carrier `q`, not a new one**: `S_eff(q)` is *the marginal of WG's own random variable*. So Q1 is a
**YES** that is entirely contained in the existing corpus. The Feshbach framing changes nothing about
the carrier.

---

## 3 · Q2 — does it improve on WG? NO `[NG]`

### 3.1 Gaussian limit: curvature = reciprocal susceptibility (exact repackaging) `[T]`
For a Gaussian measure, the coarse marginal precision after eliminating the fibre is the Schur
complement `G″ = (P A^{-1}P^T)^{-1}` and the marginal covariance is `Var = P A^{-1}P^T`; hence
`G″ · Var = I` **identically** (§7[A], error `1e-15`). This is ch02 `eq:ward` (`G″·Var=1`, Gaussian
*marginal* ⇒ identity). **In the Gaussian case the Feshbach curvature is exactly the reciprocal
susceptibility** `1/χ` — i.e. the reciprocal WG estimator. There is no independent content.

### 3.2 Off-Gaussian: `G″ ≠ 1/χ` at finite `V`, but the split is CLT noise `[T]/[V]`
The Feshbach curvature and `1/χ` are **different functionals** of the same marginal off-Gaussian:
- **one-plaquette SU(2)** (exact, `V=1`): density `ρ(z)=(2/π)√(1−z²)`, `z=½tr U`; `Var=¼` so `1/χ=4`,
  but `G″=S_eff''(0)=−ρ''(0)/ρ(0)=1` (§7[D]). Factor **4** split — the maximal finite-`V` (non-Gaussian)
  difference.
- **free composite `:φ²:`** (1D): direct histogram curvature `S_eff''` vs `1/χ` gives ratio
  `1.37 → 1.11 → 1.04` at `L=16,32,64` (§7[C]) while the excess kurtosis `κ₄ = 3.92 → 1.79 → 0.84 → 0`.
The split is a **non-Gaussianity** governed by the cumulants. By the **CLT for a mixing (gapped/
summable) field**, `q_O ⇒ N(0,χ_O)`, with `κ_n(q_O) = O(V^{1−n/2})`. Hence
```
   S_eff(q) = q²/(2χ_O) + o(1)   as V → ∞   (whenever χ_O < ∞) ,                                   (3.1)
```
so **the nonlinear Feshbach effective action is Gaussian in the thermodynamic limit and `G″→1/χ`**.
The "nonlinearity" is a `1/V`-suppressed correction with **no thermodynamic content**.

### 3.3 Ceiling, not pole `[T]`
For a **composite** (continuum KL measure, the physically relevant case for `tr F²`), neither `G″` nor
`1/χ` is the lightest mass: free `:φ²:` has `1/χ ≈ 0.24–0.25` and `S_eff'' ≈ 0.26–0.33`, both `≪
(2m)²=1` (§7[C]). This is exactly ch02's `G″ vs 1/χ` (measured 2–5% apart at β=2.5) and `χ ≠ pole`.
Both functionals equal `m₁²` **only** when the KL measure is a single atom (elementary field; §7[B]
free `φ`: `G″=m²` exactly) — the trivial case. For composites they are **ceiling-type**.

### 3.4 The full nonlinear shape is still `k=0` (thermodynamic) — blind to the edge `[NG]`
The complete `S_eff(q)` is determined by the connected correlations of `O` at **zero momentum** (the
cumulants of a volume-averaged observable). By the large-deviation / Boltzmann-equivalence picture
these are **thermodynamic** (free-energy) data of the `k=0` sector. The mass gap is the location of the
**spectral edge** (`inf supp ν`), carried by the **momentum dependence** of the two-point function —
not by any single functional of the `k=0` marginal. Formally: any functional of the `k=0` marginal is
a `k=0` object; the thin-atom construction (ch02 thin-atom Prop., `T4-4`) realises, for every `ε>0`, a
spectral measure with the **same** `k=0` cumulants to any finite order but `inf supp ν = ε → 0`. So the
**full nonlinear Feshbach**, like the variance, is **blind to the atom** — the same finite-data /
convex-order no-go, now seen to apply to the *whole* `k=0` effective action, not just `χ`. A
**nonlinear** collective coordinate `Q=f(...)` is still built from volume averages and is covered.

### 3.5 The only genuine increment: all-`k` = the standard propagator `[T-standard]/[C]/[V]`
Doing the projection at each `k` yields `C(k)=Var(q(k))=G̃(k)`. This is the **full gauge-invariant
propagator**, and its `x`-decay gives the **standard gauge-invariant effective mass** (§7[G]:
`m_eff = 0.2989` vs `m=0.3`). It does add the `k≠0` slice that WG's single number `χ_O=C(0)` omits —
**but** (i) it is the established tool (`INVENT-GAUGE-INVARIANT-MASS.md`), not new; (ii) it is an
**estimate** (`[C]/[V]`), not a lower bound — the edge (branch point of the analytic `C(k)`) is exactly
the R3 circular/infinite-data target; (iii) each `C(k)` is still a *second-moment* (Stieltjes) value of
`ν`, so finite-`k` data is finite-moment data, blind to the atom by §3.4. So the multiscale version is a
**better *estimator* of the mass than `χ_O` alone**, but not a new structural tool and not a bound.

---

## 4 · The Savvidy cell — the collective mode is **not** controlled `[NG]`

The gauge-invariant collective coordinate is `q ~ B²` (it tracks the background amplitude, since
`tr F² ≈ ½B²` on a magnetic background). On the **Savvidy one-loop cell** the effective potential is
```
   V_eff(B) = ½B² + A B²( ln(B/M²) − ½ ),   A = 11 g²/(96π²) > 0 ,     V_eff''(B) = 1 + 2A ln(B/M²) + 2A .
```
`V_eff''(B) < 0` for `B < B_c = M² exp(−(1+2A)/(2A))` — a **concave (tachyonic) region** near `B=0`,
with the Savvidy minimum at `B_min = e·B_c` on the stable side (§7[E]: sign flip, negative curvature
`V''(B_c/10)<0`). So **the gauge-invariant Feshbach on the perturbative cell inherits the physical
Savvidy instability**; the projection neither removes nor controls it. (The full-measure positivity
`1/χ_O=1.7–3.0>0` is the ensemble/compact-Haar **measure** fact — WG — not the projection: the measure
averages over all backgrounds and convexifies what the cell leaves non-convex.)

**Relation to `H1-RG-DESCENT`.** The link-space descent breaks at the soft mode because `Hess V` goes
indefinite (`η_j≥1`, `γ_j→0`, `D_j>0`, `G″≤0` at β≈2.4–2.5). The gauge-invariant Feshbach **avoids that
breakdown** — but only by replacing the mode-by-mode curvature (which it cannot keep positive) with the
**scalar marginal** `S_eff(q)`, whose curvature is `1/χ_O` (the ceiling). It does not *control* the soft
mode; it *integrates over* it and returns H1's ceiling. So it does not improve on H1 either — it trades
the H1 curvature descent for the H1 ceiling.

---

## 5 · Acid tests `[T]/[V]`

1. **Free massive scalar (known value).** `O=φ`: `S_eff` Gaussian, `G″=m²=1/Var` exactly at every `L`,
   d=1 (§7[B]); acid `m→0`: `S_eff''=m²→0` (no false floor). **Passes.**
2. **Free massless (must fail H1).** `Var(q)=1/m²→∞`, `S_eff''→0`; `q_O` not normalizable (no CLT).
   The Feshbach manufactures no floor. **Passes.**
3. **2D / free composite `tr F²` (`:φ²:`).** `χ` finite and `L`-independent for `m>0`, `→∞` for `m→0`
   (§7[C]); `f₁²=0` (atomless, continuum above threshold); `G″≈1/χ` but `≠ (2m)²`. **Passes.**
4. **One-plaquette SU(2) (exact).** `Var=¼`, `G″=1`, smooth density (atomless `q`). **Passes.**

---

## 6 · Kill-criteria ledger

| criterion | verdict | why |
|---|---|---|
| **(a) reproduces WG/H1 as repackaging** | **TRIGGERED** | `G″=1/χ` exact in the Gaussian limit (§3.1); off-Gaussian split `O(V^{-1/2})→0` by CLT and the full `S_eff` is Gaussian in the thermodynamic limit (§3.2, §7[C]); the carrier is WG's `q` (§2). |
| **(b) needs a gauge choice** | **NOT triggered** | the observable carrier `q_O` is exactly gauge-invariant, no gauge fixing (§1, §7[F]). |
| **(c) cannot control the collective mode** | **TRIGGERED** | Savvidy `V_eff''<0` inherited on the cell (§4, §7[E]); full-measure positivity is a measure fact, not the projection. |

**Net:** one hygiene criterion passed (gauge-independent), the substantive criterion (improves on WG)
failed via **(a)** and **(c)**. The proposal is a **clean reformulation**, not a new carrier and not a
gap-closing tool.

---

## 7 · Witness — `research/v3/w7/y4_feshbach.py` (`py -3.12 -I`, numpy only, one process)

Quoted (abridged):
```
[A] GAUSSIAN: || G'' * Var - I ||_max = 1.05e-15 (seed0)   => G''=1/chi exactly (repackages WG)
[B] FREE phi m=0.5: Var=4.000000 1/Var=0.250000 G''=m^2=0.250000 (all L); acid m=0.01: G''=1e-4 -> 0
[C] :phi^2: m=0.5, (2m)^2=1.0:
       L   chi=Var   1/Var   S_eff"(hist)   k4     S_eff"/(1/Var)
      16   4.11773   0.24285    0.33301    3.9233     1.3712
      32   4.07268   0.24554    0.27209    1.7871     1.1081
      64   4.06721   0.24587    0.25624    0.8364     1.0422    -> S_eff"->1/Var, k4->0 (CLT)
    exact chi(Wick) L=64: 4.10884 (1/chi=0.24338 != (2m)^2=1.0);  acid m=0.02: Var=1.97e5 -> inf
[D] ONE-PLAQUETTE SU(2): Var=0.25, 1/Var=4.000000, G''=S_eff''(0)=1.000000 (factor 4 split); atomless
[E] SAVIDDY: g^2=2.5 -> B_c=1.21e-08, B_min=3.30e-08, V''(B_c/10)=-0.1337<0, V''(B_min)=+0.0580
    toy A=0.5: V''(0.05)=-0.9957 (UNSTABLE/concave) ... V''(0.8)=+1.7769 (STABLE)  => inherits tachyon
[F] 2D SU(2) L=4 one config: |Q-Q'|/|Q| = 0.0 (gauge-inv) ; link-log |A-A'|/|A| = 3.98e-01 (DRIFTS)
[G] free m=0.3: C(k=0)=chi=11.11111=1/m^2; m_eff from C(x)=0.29889 (=m): standard eff. mass
```
Readings: **[A]** exact Gaussian repackaging; **[B]** free elementary exact (single-atom) + massless
acid; **[C]** composite: `S_eff''` and `1/χ` differ at finite `L` but converge (CLT, `κ₄→0`), both
ceiling-type (`≠ (2m)²`), atomless; **[D]** maximal finite-`V` `G″ vs 1/χ` split (one plaquette),
smooth; **[E]** Savvidy cell has `V_eff''<0` near `B=0` (inherited tachyon, not controlled); **[F]**
observable carrier exactly gauge-invariant, link-log drifts; **[G]** multiscale = standard propagator /
effective mass.

---

## 8 · Status ledger

| item | statement | label |
|---|---|---|
| gauge-independence (no gauge choice) of the observable Feshbach | `q_O` invariant under all `U→U^g`; §1/§7[F] | **[T]** |
| carrier `Var(q_O)=χ_O≠V^{-1}` | non-vacuous, but = the T4/WG carrier (not new) | **[T]**+**[V]** |
| Gaussian limit curvature = reciprocal WG | `G″·Var=I` exact, `1e-15`; ch02 `eq:ward` | **[T]**+**[V]** |
| off-Gaussian `G″≠1/χ` split | one-plaquette 1 vs 4; `:φ²:` `1.37→1.04` | **[T]**+**[V]** |
| **split collapses by CLT** | `κ₄→0`, `S_eff''→1/χ`, `S_eff` Gaussian in thermo limit | **[T]**+**[V]** |
| ceiling-not-pole for composites | `G″,1/χ ≈ 0.24–0.33 ≠ (2m)²=1` | **[T]**+**[V]** |
| **nonlinear Feshbach improves on WG** | **false** — Gaussian-limit repackaging; CLT collapse | **[NG]** |
| full `k=0` `S_eff` blind to the atom/edge | `k=0` thermodynamic object; thin-atom; R3 | **[NG]** |
| multiscale `C(k)` = standard effective mass | estimate `[C]/[V]`, not a bound, not new | **[T-standard]** |
| Savvidy cell: collective mode controlled | **false** — `V_eff''<0` inherited (tachyon) | **[NG]**+**[V]** |
| free-field acid (elementary/composite, massless) | known/0; `∞` (no floor) | **[T]**+**[V]** |
| atom `f₁²≥r₀` (+ isolation `I₂`) on the carrier | unchanged irreducible input | **[O]** |
| Clay mass gap | unbuilt | **0 %** |

---

## 9 · Sharpest line

**The gauge-invariant nonlinear Feshbach / multiscale RG projection is a correct, exactly
gauge-independent (no gauge fixing, `|ΔQ|/|Q|=0`), non-vacuous (`Var=χ_O≠V^{-1}`) reformulation — but it
is not a new carrier (it is `T4`/`WG`'s `q`), and it does not improve on WG: in the Gaussian limit its
curvature is exactly the reciprocal susceptibility `1/χ` (`G″·Var=I`, `1e-15`), and off-Gaussian the
"nonlinear" curvature `G″=S_eff''` differs from `1/χ` (one-plaquette 1 vs 4; free `:φ²:` ratio
`1.37→1.04`) only by a non-Gaussianity that **vanishes by CLT** (`κ₄→0`), so the full `k=0` effective
action is **Gaussian in the thermodynamic limit** and carries no content beyond `χ_O`; it is
ceiling-type for composites (`≠ (2m)²`), atomless, and blind to the spectral edge (thin-atom). The only
genuine increment — the all-`k` propagator `C(k)` — is the standard gauge-invariant effective mass (an
estimate, not a bound). On the Savvidy cell the gauge-invariant collective coordinate inherits the
physical tachyon (`V_eff''<0` near `B=0`); the projection does not control it. So kill criteria
**(a) repackages WG/H1** and **(c) cannot control the collective mode** are **TRIGGERED**; **(b) needs
a gauge choice** is not. The atom `f₁²≥r₀` (+ isolation) remains `[O]`. **No QED. Clay 0%.**

---

*Refs.* **Descent / Feshbach algebra:** `H1-RG-DESCENT.md` (`G″_j=(P A_j^{-1}P^T)^{-1}`, `G″=E[Σ]−D`,
`η_j` breakdown at the Savvidy soft mode β≈2.4–2.5), `S3-UV-FESHBACH.md` (§1 Schur identity, sign of `D`),
`docs/paper-v2/ch05-theorems.tex` (Feshbach–Bogoliubov lemma, Lemma NS). **Carrier / squeeze / no-go:**
`T4-CARRIER.md` (T4-1 normalization `Var(q)=χ_O` vs `Var(q̄)=χ_O/V`; T4-4 thin-atom no-go),
`WG-INVARIANT-ESTIMATOR.md` (`χ_O=0.328,0.377,0.411`), `WG-L12-L16.md` (`0.469,0.587`, slope `0.37±0.07`),
`WG-FSS-FIT.md` (saturation `[O]`), `docs/paper-v2/ch02-squeeze.tex` (`eq:chi`, `eq:squeeze`,
`eq:ward`=`G″·Var=1` Gaussian-only, §`G″ vs 1/χ` with the measured 2–5% gap, thin-atom Prop., `§schur`),
`R3-ATOM-FUNCTIONAL.md` (branch-point / infinite-data edge), `R5-EXACT-ANCHORS.md` (`:φ²:` atomless).
**Effective mass:** `INVENT-GAUGE-INVARIANT-MASS.md` (gauge-invariant spectral edge of `tr F²`, free
Maxwell `m_GI=0`). **Standard (textbook, not re-fetched):** central limit theorem for mixing
sequences / cumulant scaling `κ_n=O(V^{1−n/2})`; Savvidy one-loop effective potential (Savvidy 1977;
SU(2) `A=11g²/96π²`, negative-mode/tachyon); lattice free propagator; Källén–Lehmann. **Witness**
`research/v3/w7/y4_feshbach.py` (numpy only, one process, no downloads). No `[C]` upgraded to `[T]`;
no QED; **Clay 0%.**
